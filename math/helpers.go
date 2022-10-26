package math

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math"
	"math/big"
)

// Helpers for specialized polynomial constructs

type IntVector struct {
	Vector
}

type PolyVector struct {
	Vector
}

// ToPoly converts a coefficient vector into a polynomial and returns it.
// Warning: The coefficients must fit into an uint64!
func (v IntVector) ToPoly(baseRing *ring.Ring, isNTT bool) Polynomial {
	p := NewZeroPolynomial(baseRing)
	for i := 0; i < v.Length(); i++ {
		c := v.Element(i).(*ModInt).Value.Uint64()
		p.SetCoefficient(i, c)
	}
	p.Ref.IsNTT = isNTT
	return p
}

// Max returns the maximum integer in the vector.
func (v IntVector) Max() int64 {
	maxElement := v.Element(0).(*ModInt).Value.Int64()
	v.ForEach(func(el RingElement, _ int) {
		maxElement = max(maxElement, el.(*ModInt).Value.Int64())
	})
	return maxElement
}

// Min returns the minimum integer in the vector.
func (v IntVector) Min() int64 {
	maxElement := v.Element(0).(*ModInt).Value.Int64()
	v.ForEach(func(el RingElement, _ int) {
		maxElement = min(maxElement, el.(*ModInt).Value.Int64())
	})
	return maxElement
}

// L2Norm returns the L2 norm of an integer vector.
func (v IntVector) L2Norm() float64 {
	return math.Sqrt(float64(v.Copy().AsVec().Dot(v.Vector).(*ModInt).Uint64()))
}

// L2Norm returns the L2 norm of a polynomial vector.
func (v PolyVector) L2Norm(q *big.Int) float64 {
	acc := 0.0
	for i := 0; i < v.Length(); i++ {
		coeffs := v.Element(i).(Polynomial).Coeffs(q)
		acc += float64(coeffs.Dot(coeffs.AsVec()).(*ModInt).Value.Uint64())
	}
	return math.Sqrt(acc)
}

// InfNorm returns the infinity norm of a polynomial vector.
func (v PolyVector) InfNorm(q *big.Int) int64 {
	infNorm := v.Element(0).(Polynomial).Coeffs(q).Max()
	for i := 1; i < v.Length(); i++ {
		maxCoeff := v.Element(i).(Polynomial).Coeffs(q).Max()
		infNorm = max(maxCoeff, infNorm)
	}
	return infNorm
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
