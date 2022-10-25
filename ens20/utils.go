package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomPolynomial returns a random polynomial sampled from the given `sampler`.
func NewRandomPolynomial(baseRing *ring.Ring, sampler PolySampler) math.Polynomial {
	g := math.NewZeroPolynomial(baseRing)
	sampler.Read(g.Ref)
	return g
}

// NewRandomPolynomialVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolynomialVector(dim int, baseRing *ring.Ring, sampler PolySampler) math.Vector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		return NewRandomPolynomial(baseRing, sampler)
	})
	return v
}

// NewRandomPolynomialMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolynomialMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		return NewRandomPolynomial(baseRing, sampler)
	})
	return A
}

// NewRandomIntegerVector constructs a random vector of integers.
func NewRandomIntegerVector(dim int, mod *big.Int) math.IntVector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		return math.NewModInt(ring.RandInt(mod).Int64(), mod)
	})
	return v.AsIntVec()
}

// NewRandomTernaryIntegerVector constructs a random vector of integers where each element \in {-1, 0, 1} mod q.
func NewRandomTernaryIntegerVector(dim int, mod *big.Int) math.IntVector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		return math.NewModInt(ring.RandInt(big.NewInt(3)).Int64()-1, mod)
	})
	return v.AsIntVec()
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers.
func NewRandomIntegerMatrix(rows int, cols int, mod *big.Int) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		return math.NewModInt(ring.RandInt(mod).Int64(), mod)
	})
	return A
}

// SplitInvNTT returns the x that satisfies lhs = NTT(x_1) || NTT(x_2) || ... || NTT(x_numSplits)
func SplitInvNTT(lhs math.IntVector, numSplits int, baseRing *ring.Ring) math.PolyVector {
	lhsLen := lhs.Length()
	splitLen := lhsLen / numSplits
	return math.NewVectorFromSize(numSplits).Populate(
		func(i int) math.RingElement {
			return lhs.Slice(i*splitLen, (i+1)*splitLen).Copy().AsVec().AsIntVec().ToPoly(baseRing).InvNTT()
		}).AsPolyVec()
}

// Lmu computes the value of the function Lmu(L) = 1/k * X^mu * TrL in-place.
func Lmu(mu int, invk uint64, TrL math.Polynomial) math.Polynomial {
	// Compute X^mu
	xmu := math.NewZeroPolynomial(TrL.BaseRing).SetCoefficient(mu, 1)
	return TrL.Mul(xmu).Scale(invk).(math.Polynomial)
}

// LmuSum computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * [ \sum_{v=0}^{k-1} sig^v (f(mu, v)) ]
func LmuSum(k int, invk uint64, sig math.Automorphism, f func(int, int) math.Polynomial) math.Polynomial {
	return math.NewVectorFromSize(k).Populate(
		func(mu int) math.RingElement {
			tmp := sig.Trace(
				func(v int) math.Polynomial {
					return f(mu, v)
				}, k)
			return Lmu(mu, invk, tmp)
		}).Sum().(math.Polynomial)
}

// LmuSumOuter computes the value of the function \sum_{mu=0}^{k-1} (1/k) * X^mu * \sum_{v=0}^{k-1} \sum_{j=0}^{numSplits-1} sig^v (f(mu, v, j))
func LmuSumOuter(k, numSplits int, invk uint64, sig math.Automorphism, f func(int, int, int) math.Polynomial) math.Polynomial {
	return math.NewVectorFromSize(k).Populate(
		func(mu int) math.RingElement {
			tmp := math.NewMatrixFromDimensions(k, numSplits).Populate(
				func(v, j int) math.RingElement {
					return sig.Permute(int64(v), f(mu, v, j))
				}).Sum().(math.Polynomial)
			return Lmu(mu, invk, tmp)
		}).Sum().(math.Polynomial)
}

// CommitmentSum computes \sum_{i=0}^{k-1} \sum_{j=0}^{numSplits} alpha_{i*numSplits+j} sig^{-i} (f(i, j))
func CommitmentSum(k, numSplits int, alpha math.PolyVector, sig math.Automorphism, f func(int, int) math.Polynomial) math.Polynomial {
	return math.NewMatrixFromDimensions(k, numSplits).Populate(
		func(i int, j int) math.RingElement {
			index := (i*numSplits + j) % alpha.Length()
			return sig.Permute(int64(-i), f(i, j)).Mul(alpha.Element(index))
		}).Sum().(math.Polynomial)
}
