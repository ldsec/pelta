package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
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
	return logging.LogShortExecution("LmuSumOuter", "calculating", func() interface{} {
		out := fastmath.NewPolyVec(k, params.config.BaseRing).NTT()
		for mu := 0; mu < k; mu++ {
			func(k, mu int) {
				tmp2 := fastmath.NewPoly(params.config.BaseRing)
				for v := 0; v < k; v++ {
					gen := params.config.ValueCache.Get("exp", int64(v), func() uint64 { return params.Sig.Exponent(int64(v)) })
					for j := 0; j < numSplits; j++ {
						sigf := f(mu, v, j).InvNTT().PermuteWithGen(gen)
						tmp2.Add(sigf)
					}
				}
				out.Set(mu, Lmu(mu, invk, tmp2.NTT(), params))
			}(k, mu)
		}
		return out.Sum()
	}).(*fastmath.PolyNTT)
	// 	tmp := fastmath.NewPolyVec(k, params.config.BaseRing).NTT()
	// 	tmp.Populate(func(mu int) *fastmath.PolyNTT {
	// 		tmp2 := fastmath.NewPolyMatrix(k, numSplits, params.config.BaseRing).NTT()
	// 		tmp2.Populate(func(v, j int) *fastmath.PolyNTT {
	// 			return f(mu, v, j)
	// 		})
	// 		tmp2.Update(func(v, j int, old *fastmath.PolyNTT) *fastmath.PolyNTT {
	// 			gen := params.config.ValueCache.Get("exp", int64(v), func() uint64 { return params.Sig.Exponent(int64(v)) })
	// 			return old.InvNTT().PermuteWithGen(gen).NTT()
	// 		})
	// 		innerSum := tmp2.Sum()
	// 		return Lmu(mu, invk, innerSum, params)
	// 	})
	// 	return tmp.Sum()
	// }).(*fastmath.PolyNTT)
}

// CommitmentSum computes \sum_{i=0}^{k-1} \sum_{j=0}^{numSplits} alpha_{i*numSplits+j} sig^{-i} (f(i, j))
func CommitmentSum(k, numSplits int, alpha *fastmath.PolyNTTVec, f func(int, int) *fastmath.PolyNTT, params PublicParams) *fastmath.PolyNTT {
	return logging.LogShortExecution("CommitmentSum", "calculating", func() interface{} {
		presum := fastmath.NewPolyMatrix(k, numSplits, params.config.BaseRing).NTT()
		presum.Populate(func(i, j int) *fastmath.PolyNTT {
			fOut := f(i, j)
			gen := params.config.ValueCache.Get("exp", int64(-i), func() uint64 { return params.Sig.Exponent(int64(-i)) })
			return fOut.InvNTT().PermuteWithGen(gen).NTT()
		})
		presum.Update(func(i, j int, old *fastmath.PolyNTT) *fastmath.PolyNTT {
			index := (i*numSplits + j) % alpha.Size()
			return old.Mul(alpha.Get(index))
		})
		return presum.Sum()
	}).(*fastmath.PolyNTT)
}
