package math

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math"
	"math/big"
)

// Helpers for specialized polynomial constructs

type PolyArray struct {
	*MultiArray
}

type IntVector struct {
	Vector
}

type PolyVector struct {
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

// L2Norm returns the L2 norm of an integer vector.
func (v IntVector) L2Norm() float64 {
	return math.Sqrt(float64(v.Dot(v.Vector).(*ModInt).Uint64()))
}

// L2Norm returns the L2 norm of a polynomial vector.
func (v PolyVector) L2Norm(q *big.Int) float64 {
	acc := 0.0
	for i := 0; i < v.Length(); i++ {
		acc += v.Element(i).(Polynomial).Coeffs(q).L2Norm()
	}
	return math.Sqrt(acc)
}

// -- Conversion helpers

// AsMatrix converts the representation to a matrix.
func (m *MultiArray) AsMatrix() Matrix {
	// assert len(a.Dimensions()) == 2
	if len(m.Dimensions()) != 2 {
		panic(fmt.Sprintf("AsMatrix: Cannot convert a multi-array with %d dimensions into a matrix", m.Dimensions()))
	}
	return Matrix{m}
}

// AsVec converts the representation to a vector.
func (m *MultiArray) AsVec() Vector {
	// assert len(a.Dimensions()) == 1
	if len(m.Dimensions()) != 1 {
		panic(fmt.Sprintf("AsVec: Cannot convert a multi-array with %d dimensions into a vector", m.Dimensions()))
	}
	return Vector{m}
}

// AsPolyArray converts to the representation of a multi array of polynomials.
// Allows calling specialized methods.
func (m *MultiArray) AsPolyArray() PolyArray {
	// assert type
	if _, ok := m.ElementAtIndex(0).(Polynomial); !ok {
		panic(fmt.Sprintf("AsPolyArray: Not a polynomial array"))
	}
	return PolyArray{m}
}

// AsIntVec converts to the representation of a mod int.
// Allows calling specialized methods.
func (v Vector) AsIntVec() IntVector {
	// assert type
	if _, ok := v.Element(0).(*ModInt); !ok {
		panic(fmt.Sprintf("AsIntVec: Not an integer vector"))
	}
	return IntVector{v}
}

// AsPolyVec converts to the representation of a vector of polynomials.
// Allows calling specialized methods.
func (v Vector) AsPolyVec() PolyVector {
	if _, ok := v.Element(0).(Polynomial); !ok {
		panic(fmt.Sprintf("AsPolyVec: Not a polynomial vector"))
	}
	return PolyVector{v}
}
