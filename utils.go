package main

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/lattigo/v2/ring"
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
		return math.NewModInt(ring.RandInt(mod), mod)
	})
	return v
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers.
func NewRandomIntegerMatrix(rows int, cols int, mod *big.Int) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		return math.NewModInt(ring.RandInt(mod), mod)
	})
	return A
}

func Sig(exp int, p math.Polynomial) math.Polynomial {
	// TODO complete
	return p
}
