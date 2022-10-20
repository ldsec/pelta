package main

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

// NewRandomIntegerVector constructs a random 2D vector of integers.
func NewRandomIntegerVector(dim int, mod *big.Int) math.Vector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		return math.NewModInt(ring.RandInt(mod).Int64(), mod)
	})
	return v
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers.
func NewRandomIntegerMatrix(rows int, cols int, mod *big.Int) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		return math.NewModInt(ring.RandInt(mod).Int64(), mod)
	})
	return A
}

// SplitInvNTT returns the x that satisfies lhs = NTT(x_1) || NTT(x_2) || ... || NTT(x_numSplits)
func SplitInvNTT(lhs math.Vector, numSplits, d int, baseRing *ring.Ring) math.Vector {
	return math.NewVectorFromSize(numSplits).Populate(
		func(i int) math.RingElement {
			return lhs.Slice(i*d, (i+1)*d).Copy().AsVec().AsCoeffs().ToPoly(baseRing).InvNTT()
		})
}

// Lmu computes the value of the function Lmu(L) = 1/k * X^mu * TrL
func Lmu(mu int, TrL math.Polynomial, invk *math.ModInt) math.Polynomial {
	return TrL.Copy().(math.Polynomial).
		RRot(mu).
		Scale(invk.Uint64()).(math.Polynomial)
}
