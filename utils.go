package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// RandomVector constructs a vector, whose elements sampled from the given `sampler`.
func RandomVector(dim int, baseRing *ring.Ring, sampler PolySampler) MultiArray {
	v := NewMultiArray([]int{dim}, baseRing)
	for i := 0; i < dim; i++ {
		sampler.Read(v.ElementAtIndex(i))
	}
	return v
}

// RandomMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func RandomMatrix(rows int, cols int, baseRing *ring.Ring, sampler PolySampler) MultiArray {
	A := NewMultiArray([]int{rows, cols}, baseRing)
	for i := 0; i < rows*cols; i++ {
		sampler.Read(A.ElementAtIndex(i))
	}
	return A
}

// DotProduct performs a dot product between two vectors.
func DotProduct(a MultiArray, b MultiArray, multOperator BinaryOperator, baseRing *ring.Ring) *ring.Poly {
	if !a.IsVector() || !b.IsVector() {
		panic(fmt.Errorf("trying to take the dot product among non-vectors"))
	}
	c := baseRing.NewPoly()
	for i := 0; i < a.Length(); i++ {
		tmp := multOperator.Apply(a.ElementAtIndex(i), b.ElementAtIndex(i), baseRing)
		baseRing.Add(c, tmp, c)
	}
	return c
}

// MatrixVectorMul performs a matrix vector multiplication and returns the result Ax = y
func MatrixVectorMul(A MultiArray, x MultiArray, multOperator BinaryOperator, baseRing *ring.Ring) []*ring.Poly {
	y := make([]*ring.Poly)
	for i := 0; i < len(A); i++ {
		y[i] = DotProduct(A[i], x, multOperator, baseRing)
	}
	return y
}

// ApplyMatrix applies a given operator to the given matrix.
func ApplyMatrix(M [][]*ring.Poly, op UnaryOperator, baseRing *ring.Ring) [][]*ring.Poly {
	out := make([][]*ring.Poly, len(M))
	for i := 0; i < len(out); i++ {
		out[i] = make([]*ring.Poly, len(M[0]))
		for j := 0; j < len(out[i]); j++ {
			out[i][j] = op.Apply(M[i][j], baseRing)
		}
	}
	return out
}

// VectorAdd adds the two given vectors, whose elements are added coefficient wise.
func VectorAdd(a []*ring.Poly, b []*ring.Poly, baseRing *ring.Ring) []*ring.Poly {
	c := make([]*ring.Poly, len(a))
	for i := 0; i < len(c); i++ {
		c[i] = baseRing.NewPoly()
		baseRing.Add(a[i], b[i], c[i])
	}
	return c
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
