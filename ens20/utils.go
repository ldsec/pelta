package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// SplitInvNTT returns the x that satisfies lhs = NTT(x_1) || NTT(x_2) || ... || NTT(x_numSplits)
func SplitInvNTT(lhs rings.ZIntVector, numSplits int, mod *big.Int, baseRing *ring.Ring) rings.PolyVector {
	lhsLen := lhs.Length()
	splitLen := lhsLen / numSplits
	return rings.NewPolyVec(algebra.NewVectorFromSize(numSplits).Populate(
		func(i int) algebra.Element {
			return rings.NewZIntVec(lhs.Slice(i*splitLen, (i+1)*splitLen).Copy().AsVec()).ToPoly(baseRing, mod, true).InvNTT()
		}))
}

// Lmu computes the value of the function Lmu(L) = 1/k * X^mu * TrL in-place.
func Lmu(mu int, invk uint64, TrL rings.Polynomial) rings.Polynomial {
	// Compute X^mu
	xmu := rings.NewZeroPolynomial(TrL.BaseRing).SetCoefficient(mu, 1)
	return TrL.Mul(xmu).Scale(invk).(rings.Polynomial)
}

// LmuSum computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * [ \sum_{v=0}^{k-1} sig^v (f(mu, v)) ]
func LmuSum(k int, invk uint64, sig math.Automorphism, f func(int, int) rings.Polynomial) rings.Polynomial {
	return algebra.NewVectorFromSize(k).Populate(
		func(mu int) algebra.Element {
			tmp := sig.Trace(
				func(v int) rings.Polynomial {
					return f(mu, v)
				}, k)
			return Lmu(mu, invk, tmp)
		}).Sum().(rings.Polynomial)
}

// LmuSumOuter computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * \sum_{v=0}^{k-1} \sum_{j=0}^{numSplits-1} sig^v (f(mu, v, j))
func LmuSumOuter(k, numSplits int, invk uint64, sig math.Automorphism, f func(int, int, int) rings.Polynomial) rings.Polynomial {
	return algebra.NewVectorFromSize(k).Populate(
		func(mu int) algebra.Element {
			tmp := algebra.NewMatrixFromDimensions(k, numSplits).Populate(
				func(v, j int) algebra.Element {
					return sig.Permute(int64(v), f(mu, v, j))
				}).Sum().(rings.Polynomial)
			return Lmu(mu, invk, tmp)
		}).Sum().(rings.Polynomial)
}

// CommitmentSum computes \sum_{i=0}^{k-1} \sum_{j=0}^{numSplits} alpha_{i*numSplits+j} sig^{-i} (f(i, j))
func CommitmentSum(k, numSplits int, alpha rings.PolyVector, sig math.Automorphism, f func(int, int) rings.Polynomial) rings.Polynomial {
	unscaled := algebra.NewMatrixFromDimensions(k, numSplits).Populate(
		func(i int, j int) algebra.Element {
			return sig.Permute(int64(-i), f(i, j))
		})
	return unscaled.Map(
		func(el algebra.Element, i int, j int) algebra.Element {
			index := (i*numSplits + j) % alpha.Length()
			return el.Mul(alpha.Element(index))
		}).Sum().(rings.Polynomial)
}
