package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomVector(dim int, baseRing *ring.Ring, sampler PolySampler) MultiArray {
	v := NewMultiArray([]int{dim}, baseRing)
	for i := 0; i < dim; i++ {
		sampler.Read(v.ElementAtIndex(i))
	}
	return v
}

// NewRandomMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) MultiArray {
	A := NewMultiArray([]int{cols, rows}, baseRing)
	for i := 0; i < rows*cols; i++ {
		sampler.Read(A.ElementAtIndex(i))
	}
	return A
}

// DotProduct performs a dot product between two vectors.
func DotProduct(a MultiArray, b MultiArray, multOperator BinaryOperator, baseRing *ring.Ring) *ring.Poly {
	if !a.IsVector() || !b.IsVector() {
		panic(fmt.Errorf("invalid product"))
	}
	c := baseRing.NewPoly()
	for i := 0; i < a.Length(); i++ {
		// TODO: convert to muladd
		tmp := baseRing.NewPoly()
		multOperator.Apply(a.ElementAtIndex(i), b.ElementAtIndex(i), tmp, baseRing)
		baseRing.Add(c, tmp, c)
	}
	return c
}

// MatrixVectorMul performs a matrix vector multiplication and returns the result Ax = y
func MatrixVectorMul(A MultiArray, x MultiArray, multOperator BinaryOperator, baseRing *ring.Ring) MultiArray {
	if !A.IsMatrix() || !x.IsVector() {
		panic(fmt.Errorf("invalid product"))
	}
	outputDim := A.Dimensions()[0]
	y := NewMultiArray([]int{outputDim}, baseRing)
	for i := 0; i < A.Dimensions()[1]; i++ {
		y.SetElementAtIndex(i, DotProduct(A.MatrixRowSlice(i), x, multOperator, baseRing))
	}
	return y
}

// VectorAdd adds the two given vectors, whose elements are added coefficient wise.
func VectorAdd(a MultiArray, b MultiArray, baseRing *ring.Ring) MultiArray {
	if !a.IsVector() || !b.IsVector() {
		panic(fmt.Errorf("invalid sum"))
	}
	c := make([]*ring.Poly, len(a.Array))
	for i := 0; i < len(c); i++ {
		c[i] = baseRing.NewPoly()
		baseRing.Add(a.ElementAtIndex(i), b.ElementAtIndex(i), c[i])
	}
	return NewVectorFromSlice(c)
}

// VectorSum sums the polynomials in the given vector.
func VectorSum(v []*ring.Poly, baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	for i := 0; i < len(v); i++ {
		baseRing.Add(out, v[i], out)
	}
	return out
}

// MatrixTranspose transposes a given 2x2 matrix.
func MatrixTranspose(Ain [][]*ring.Poly) [][]*ring.Poly {
	lenX := len(Ain[0])
	lenY := len(Ain)
	Aout := make([][]*ring.Poly, lenX)
	for i := 0; i < lenX; i++ {
		Aout[i] = make([]*ring.Poly, lenY)
	}
	for i := 0; i < lenX; i++ {
		for j := 0; j < lenY; j++ {
			Aout[i][j] = Ain[j][i]
		}
	}
	return Aout
}
