package rings

import (
	"fmt"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// Helpers for specialized polynomial constructs

type ZIntVector struct {
	algebra.Vector
}

type PolyVector struct {
	algebra.Vector
}

// ToPoly converts a coefficient vector into a polynomial and returns it.
// Warning: The coefficients must fit into an uint64!
func (v ZIntVector) ToPoly(baseRing *ring.Ring, mod *big.Int, isNTT bool) Polynomial {
	p := NewZeroPolynomial(baseRing)
	for i := 0; i < v.Length(); i++ {
		c := v.Element(i).(ZInt).Uint64WithMod(mod)
		p.SetCoefficient(i, c)
	}
	p.Ref.IsNTT = isNTT
	return p
}

// Max returns the maximum integer in the vector.
func (v ZIntVector) Max() int64 {
	maxElement := v.Element(0).(ZInt).Int64()
	v.ForEach(func(el algebra.Element, _ int) {
		maxElement = max(maxElement, el.(ZInt).Int64())
	})
	return maxElement
}

// Min returns the minimum integer in the vector.
func (v ZIntVector) Min() int64 {
	maxElement := v.Element(0).(ZInt).Int64()
	v.ForEach(func(el algebra.Element, _ int) {
		maxElement = min(maxElement, el.(ZInt).Int64())
	})
	return maxElement
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

// NewZIntVec converts to the representation of an int vector.
// Allows calling specialized methods.
func NewZIntVec(v algebra.Vector) ZIntVector {
	// assert type
	if _, ok := v.Element(0).(ZInt); !ok {
		panic(fmt.Sprintf("NewZIntVec: Not an integer vector"))
	}
	return ZIntVector{v}
}

// NewPolyVec converts to the representation of a vector of polynomials.
// Allows calling specialized methods.
func NewPolyVec(v algebra.Vector) PolyVector {
	if _, ok := v.Element(0).(Polynomial); !ok {
		panic(fmt.Sprintf("NewPolyVec: Not a polynomial vector"))
	}
	return PolyVector{v}
}
