package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// SplitInvNTT extracts the underlying polynomials from an integer vector such that the integer vector
// satisfies lhs = NTT(x_1) || NTT(x_2) || ... || NTT(x_numSplits)
func SplitInvNTT(lhs fastmath.IntVec, baseRing *ring.Ring) fastmath.PolyVec {
	splits := fastmath.NewPolyVec(len(lhs.UnderlyingPolys()), baseRing)
	splits.Populate(func(i int) fastmath.Poly {
		split := lhs.UnderlyingPolys()[i].Copy()
		return *split.InvNTT()
	})
	return splits
}

// Lmu computes the value of the function Lmu(L) = 1/k * X^mu * TrL in-place.
func Lmu(mu int, invk uint64, TrL *fastmath.Poly, baseRing *ring.Ring) fastmath.Poly {
	// Compute X^mu
	xmu := fastmath.NewZeroPoly(baseRing)
	xmu.Set(mu, 1)
	xmu.NTT()
	return *TrL.MulCoeffs(&xmu).Scale(invk)
}

// LmuSum computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * [ \sum_{v=0}^{k-1} sig^v (f(mu, v)) ]
func LmuSum(k int, invk uint64, sig fastmath.Automorphism, f func(int, int) fastmath.Poly, baseRing *ring.Ring) fastmath.Poly {
	tmp := fastmath.NewPolyVec(k, baseRing)
	tmp.Populate(func(mu int) fastmath.Poly {
		traceFunc := func(v int) fastmath.Poly {
			return f(mu, v)
		}
		traceResult := fastmath.Trace(sig, traceFunc, k, baseRing)
		return Lmu(mu, invk, traceResult.NTT(), baseRing)
	})
	return tmp.Sum()
}

// LmuSumOuter computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * \sum_{v=0}^{k-1} \sum_{j=0}^{numSplits-1} sig^v (f(mu, v, j))
func LmuSumOuter(k, numSplits int, invk uint64, sig fastmath.Automorphism, f func(int, int, int) fastmath.Poly, baseRing *ring.Ring) fastmath.Poly {
	tmp := fastmath.NewPolyVec(k, baseRing)
	tmp.Populate(func(mu int) fastmath.Poly {
		tmp2 := fastmath.NewPolyMatrix(k, numSplits, baseRing)
		tmp2.Populate(func(v, j int) fastmath.Poly {
			fOut := f(mu, v, j)
			return fOut.Permute(int64(v), sig)
		})
		innerSum := tmp2.Sum()
		return Lmu(mu, invk, innerSum.NTT(), baseRing)
	})
	return tmp.Sum()
}

// CommitmentSum computes \sum_{i=0}^{k-1} \sum_{j=0}^{numSplits} alpha_{i*numSplits+j} sig^{-i} (f(i, j))
func CommitmentSum(k, numSplits int, alpha fastmath.PolyVec, sig fastmath.Automorphism, f func(int, int) fastmath.Poly, baseRing *ring.Ring) fastmath.Poly {
	presum := fastmath.NewPolyMatrix(k, numSplits, baseRing)
	presum.Populate(func(i, j int) fastmath.Poly {
		fOut := f(i, j)
		return fOut.Permute(int64(-i), sig)
	})
	presum.Update(func(i, j int, old fastmath.Poly) fastmath.Poly {
		index := (i*numSplits + j) % alpha.Size()
		return *old.NTT().MulCoeffs(alpha.Get(index))
	})
	return presum.Sum()
}
