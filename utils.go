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
func NewRandomVector(dim int, baseRing *ring.Ring, sampler PolySampler) math.MultiArray {
	v := math.NewMultiArray([]int{dim}, baseRing)
	for i := 0; i < dim; i++ {
		sampler.Read(v.ElementAtIndex(i))
	}
	return v
}

// NewRandomMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) math.MultiArray {
	A := math.NewMultiArray([]int{cols, rows}, baseRing)
	for i := 0; i < rows*cols; i++ {
		sampler.Read(A.ElementAtIndex(i))
	}
	return A
}
