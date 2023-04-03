package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"sync"
)

// SplitInvNTT extracts the underlying polynomials from an integer vector such that the integer vector
// satisfies lhs = NTT(x_1) || NTT(x_2) || ... || NTT(x_k)
func SplitInvNTT(lhs *fastmath.IntVec, params PublicParams) *fastmath.PolyVec {
	splits := fastmath.NewPolyVec(len(lhs.UnderlyingPolys()), params.config.BaseRing)
	splits.Populate(func(i int) *fastmath.Poly {
		// We assume that LHS is in NTT form.
		split := fastmath.ForceNTT(lhs.UnderlyingPolys()[i].Copy())
		// Move back to the poly space.
		return split.InvNTT()
	})
	return splits
}

// Lmu computes the value of the function Lmu(L) = 1/k * X^mu * TrL in-place.
func Lmu(mu int, invk uint64, Tr *fastmath.PolyNTT, params PublicParams) *fastmath.PolyNTT {
	// Compute X^mu
	xmu := fastmath.NewPoly(params.config.BaseRing)
	xmu.SetForce(mu, 1)
	return Tr.Mul(xmu.NTT()).Scale(invk)
}

// LmuSum computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * [ \sum_{v=0}^{k-1} sig^v (f(mu, v)) ]
func LmuSum(k int, invk uint64, f func(int, int) *fastmath.PolyNTT, params PublicParams) *fastmath.PolyNTT {
	return logging.LogShortExecution("LmuSum", "calculating", func() interface{} {
		tmp := fastmath.NewPolyVec(k, params.config.BaseRing).NTT()
		tmp.Populate(func(mu int) *fastmath.PolyNTT {
			traceFunc := func(v int) *fastmath.Poly {
				fOut := f(mu, v)
				return fOut.InvNTT()
			}
			traceResult := fastmath.Trace(params.Sig, traceFunc, k, params.config.BaseRing)
			return Lmu(mu, invk, traceResult.NTT(), params)
		})
		return tmp.Sum()
	}).(*fastmath.PolyNTT)
}

// LmuSumOuter computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * \sum_{v=0}^{k-1} \sum_{j=0}^{numSplits-1} sig^v (f(mu, v, j))
func LmuSumOuter(k, numSplits int, invk uint64, f func(int, int, int) *fastmath.PolyNTT, params PublicParams) *fastmath.PolyNTT {
	// set the number of batches for the operation
	// should ideally correspond to the # of cores
	numBatches := 16
	return logging.LogShortExecution("LmuSumOuter", "calculating", func() interface{} {
		out := fastmath.NewPolyVec(k, params.config.BaseRing).NTT()
		for mu := 0; mu < k; mu++ {
			tmp2 := fastmath.NewPoly(params.config.BaseRing)
			for v := 0; v < k; v++ {
				gen := params.config.Cache.SigmaExpCache[int64(v)]
				sigfs := fastmath.NewPolyVec(numSplits, params.config.BaseRing)
				opPerBatch := numSplits / numBatches
				// split batches across different goroutines
				var wg sync.WaitGroup
				wg.Add(numBatches + 1)
				for b := 0; b < numBatches+1; b++ {
					go func(batchId int) {
						for j := batchId * opPerBatch; j < (batchId+1)*opPerBatch; j++ {
							// for the last batch
							if j >= numSplits {
								break
							}
							sigf := f(mu, v, j).InvNTT().PermuteWithGen(gen)
							sigfs.Set(j, sigf)
						}
						wg.Done()
					}(b)
				}
				wg.Wait()
				tmp2.Add(sigfs.Sum())
			}
			out.Set(mu, Lmu(mu, invk, tmp2.NTT(), params))
		}
		return out.Sum()
	}).(*fastmath.PolyNTT)
}

// LmuSumOuter computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * \sum_{v=0}^{k-1} \sum_{j=0}^{numSplits-1} sig^v (f(mu, v, j))
func LmuSumOuterDot(k, numSplits int, invk uint64, f1 func(int, int, int) *fastmath.PolyNTTVec, f2 func(int, int, int) *fastmath.PolyNTTVec, params PublicParams) *fastmath.PolyNTT {
	return logging.LogShortExecution("LmuSumOuterDot", "calculating", func() interface{} {
		out := fastmath.NewPolyVec(k, params.config.BaseRing).NTT()
		for mu := 0; mu < k; mu++ {
			for v := 0; v < k; v++ {
				gen := params.config.Cache.SigmaExpCache[int64(v)]
				tmp2 := fastmath.NewPoly(params.config.BaseRing).NTT()
				for j := 0; j < numSplits; j++ {
					v1 := f1(mu, v, j)
					v1.Update(func(_ int, old *fastmath.PolyNTT) *fastmath.PolyNTT {
						return old.InvNTT().PermuteWithGen(gen).NTT()
					})
					v2 := f2(mu, v, j)
					v2.Update(func(_ int, old *fastmath.PolyNTT) *fastmath.PolyNTT {
						return old.InvNTT().PermuteWithGen(gen).NTT()
					})
					tmp2.Add(v1.Dot(v2))
				}
				out.Set(mu, Lmu(mu, invk, tmp2, params))
			}
		}
		return out.Sum()
	}).(*fastmath.PolyNTT)
}

// CommitmentSum computes \sum_{i=0}^{k-1} \sum_{j=0}^{numSplits} alpha_{i*numSplits+j} sig^{-i} (f(i, j))
func CommitmentSum(k, numSplits int, alpha *fastmath.PolyNTTVec, f func(int, int) *fastmath.PolyNTT, params PublicParams) *fastmath.PolyNTT {
	return logging.LogShortExecution("CommitmentSum", "calculating", func() interface{} {
		presum := fastmath.NewPolyMatrix(k, numSplits, params.config.BaseRing).NTT()
		presum.Populate(func(i, j int) *fastmath.PolyNTT {
			fOut := f(i, j)
			gen := params.config.Cache.SigmaExpCache[int64(-i)]
			return fOut.InvNTT().PermuteWithGen(gen).NTT()
		})
		presum.Update(func(i, j int, old *fastmath.PolyNTT) *fastmath.PolyNTT {
			index := (i*numSplits + j) % alpha.Size()
			return old.Mul(alpha.Get(index))
		})
		return presum.Sum()
	}).(*fastmath.PolyNTT)
}
