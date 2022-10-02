package main

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/lattigo/v2/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomVector(dim int, baseRing *ring.Ring, sampler PolySampler) math.Vector {
	v := math.NewVectorFromSize(dim, baseRing)
	v.ForEach(func(el math.RingElement, _ int) {
		sampler.Read(el.(math.Polynomial).Ref)
	})
	return v
}

// NewRandomMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) math.Matrix {
	A := math.NewMatrixFromDimensions(rows, cols, baseRing)
	A.ForEach(func(el math.RingElement, _ int, _ int) {
		sampler.Read(el.(math.Polynomial).Ref)
	})
	return A
}
