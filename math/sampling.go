package math

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomPolynomial returns a random polynomial sampled from the given `sampler`.
func NewRandomPolynomial(baseRing *ring.Ring, sampler PolySampler) rings.Polynomial {
	g := rings.NewZeroPolynomial(baseRing)
	sampler.Read(g.Ref)
	return g
}

// NewRandomPolynomialVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolynomialVector(dim int, baseRing *ring.Ring, sampler PolySampler) algebra.Vector {
	v := algebra.NewVectorFromSize(dim).Populate(
		func(_ int) algebra.Element {
			return NewRandomPolynomial(baseRing, sampler)
		})
	return v
}

// NewRandomPolynomialMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolynomialMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) algebra.Matrix {
	A := algebra.NewMatrixFromDimensions(rows, cols).Populate(
		func(_, _ int) algebra.Element {
			return NewRandomPolynomial(baseRing, sampler)
		})
	return A
}

// NewRandomIntegerVector constructs a random vector of integers mod n.
func NewRandomIntegerVector(dim int, n *big.Int) rings.ZIntVector {
	v := algebra.NewVectorFromSize(dim).Populate(
		func(_ int) algebra.Element {
			return rings.NewZqInt(ring.RandInt(n).Int64(), n)
		})
	return rings.NewZIntVec(v)
}

// NewRandomTernaryIntegerVector constructs a random vector of integers where each element \in {-1, 0, 1}.
func NewRandomTernaryIntegerVector(dim int, n *big.Int) rings.ZIntVector {
	v := algebra.NewVectorFromSize(dim).Populate(
		func(_ int) algebra.Element {
			return rings.NewZInt(ring.RandInt(big.NewInt(3)).Int64() - 1)
		})
	return rings.NewZIntVec(v)
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers mod n.
func NewRandomIntegerMatrix(rows int, cols int, n *big.Int) algebra.Matrix {
	A := algebra.NewMatrixFromDimensions(rows, cols).Populate(
		func(_, _ int) algebra.Element {
			return rings.NewZqInt(ring.RandInt(n).Int64(), n)
		})
	return A
}
