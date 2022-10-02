package main

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/lattigo/v2/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomPolynomialVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolynomialVector(dim int, baseRing *ring.Ring, sampler PolySampler) math.Vector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		tmp := math.NewPolynomial(baseRing)
		sampler.Read(tmp.Ref)
		return tmp
	})
	return v
}

// NewRandomPolynomialMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolynomialMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		tmp := math.NewPolynomial(baseRing)
		sampler.Read(tmp.Ref)
		return tmp
	})
	return A
}

// NewRandomIntegerVector constructs a random 2D vector of integers.
func NewRandomIntegerVector(dim int, mod int) math.Vector {
	v := math.NewVectorFromSize(dim).Populate(func(_ int) math.RingElement {
		// TODO fix
		return math.NewModInt(1, mod)
	})
	return v
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers.
func NewRandomIntegerMatrix(rows int, cols int, mod int) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols).Populate(func(_, _ int) math.RingElement {
		// TODO fix
		return math.NewModInt(1, mod)
	})
	return A
}
