package math

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"math"
)

// Helpers for specialized polynomial constructs

type PolyArray struct {
	*MultiArray
}

type IntVector struct {
	Vector
}

// NTT converts the vector of polynomials into the NTT space in-place.
// p_i => NTT(p_i)
func (m PolyArray) NTT() PolyArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).NTT()
	})
	return m
}

// InvNTT converts the vector of NTT polynomials back into the poly space in-place.
// p_i => InvNTT(p_i)
func (m PolyArray) InvNTT() PolyArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).InvNTT()
	})
	return m
}

func (m PolyArray) MulPoly(p Polynomial) PolyArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).Mul(p)
	})
	return m
}

// ToPoly converts a coefficient vector into a polynomial.
// Warning: The coefficients must fit into an uint64!
func (v IntVector) ToPoly(baseRing *ring.Ring) Polynomial {
	p := NewZeroPolynomial(baseRing)
	for i := 0; i < v.Length(); i++ {
		c := v.Element(i).(*ModInt).Value.Uint64()
		p.SetCoefficient(i, c)
	}
	return p
}

// L2Norm returns the L2 norm of the vector.
func (v IntVector) L2Norm() float64 {
	return math.Sqrt(float64(v.Dot(v.Vector).(*ModInt).Uint64()))
}

// -- Conversion helpers

// AsMatrix converts the representation to a matrix.
func (m *MultiArray) AsMatrix() Matrix {
	// assert len(a.Dimensions()) == 2
	return Matrix{m}
}

// AsVec converts the representation to a vector.
func (m *MultiArray) AsVec() Vector {
	// assert len(a.Dimensions()) == 1
	return Vector{m}
}

// AsPolyArray converts the representation of a multi array of polynomials.
// Allows calling specialized methods.
func (m *MultiArray) AsPolyArray() PolyArray {
	// assert type
	return PolyArray{m}
}

// AsCoeffs converts a mod int vector into a coefficient vector.
func (v Vector) AsCoeffs() IntVector {
	// assert type
	return IntVector{v}
}
