package math

import "github.com/ldsec/lattigo/v2/ring"

// -- Polynomial helpers

// Add adds two polynomials in-place.
// p, q => p += q
func (p Polynomial) Add(q Polynomial, baseRing *ring.Ring) Polynomial {
	p.ApplyBinaryOp(SimpleAdd{}, q, p, baseRing)
	return p
}

// NTT converts into the NTT space in-place.
// p => p := NTT(p)
func (p Polynomial) NTT(baseRing *ring.Ring) Polynomial {
	p.ApplyUnaryOp(NTT{}, p, baseRing)
	return p
}

// InvNTT converts back into the poly space in-place.
// p => p := InvNTT(p)
func (p Polynomial) InvNTT(baseRing *ring.Ring) Polynomial {
	p.ApplyUnaryOp(InvNTT{}, p, baseRing)
	return p
}

// Scale scales the coefficients of the polynomial in-place.
// p => p := c*p
func (p Polynomial) Scale(factor uint64, baseRing *ring.Ring) Polynomial {
	p.ApplyUnaryOp(Scale{factor}, p, baseRing)
	return p
}

// Mul multiplies two polynomials coefficient-wise in-place.
// p, q => p *= q
func (p Polynomial) Mul(q Polynomial, baseRing *ring.Ring) Polynomial {
	p.ApplyBinaryOp(NTTMul{}, q, p, baseRing)
	return p
}

// Neg negates the polynomial in-place.
func (p Polynomial) Neg(baseRing *ring.Ring) Polynomial {
	// TODO implement.
	return p
}

// -- Multiarray helpers

// Add adds two arrays coefficient-wise in-place.
// p, q => p += q
func (a *MultiArray) Add(q *MultiArray, baseRing *ring.Ring) *MultiArray {
	a.ApplyBinaryOp(SimpleAdd{}, q, a, baseRing)
	return a
}

// Sum sums the polynomials in the array, returning the result.
func (a *MultiArray) Sum(baseRing *ring.Ring) Polynomial {
	out := NewPolynomial(baseRing)
	a.ApplyReduction(SimpleAdd{}, out, baseRing)
	return out
}

// Product takes the coefficient-wise product of all the elements in the array, returning the result.
func (a *MultiArray) Product(baseRing *ring.Ring) Polynomial {
	out := NewPolynomial(baseRing)
	a.ApplyReduction(NTTMul{}, out, baseRing)
	return out
}

// NTT converts the vector into the NTT space in-place.
// p => p := NTT(p)
func (a *MultiArray) NTT(baseRing *ring.Ring) *MultiArray {
	a.ApplyUnaryOp(NTT{}, a, baseRing)
	return a
}

// InvNTT converts the elements back into the poly space in-place.
// p => p := InvNTT(p)
func (a *MultiArray) InvNTT(baseRing *ring.Ring) *MultiArray {
	a.ApplyUnaryOp(InvNTT{}, a, baseRing)
	return a
}

// Scale scales the coefficients of the polynomials in-place.
// p => p := c*p
func (a *MultiArray) Scale(factor uint64, baseRing *ring.Ring) *MultiArray {
	a.ApplyUnaryOp(Scale{factor}, a, baseRing)
	return a
}

// Mul multiplies the elements coefficient-wise in-place.
// p, q => p *= q
func (a *MultiArray) Mul(q *MultiArray, baseRing *ring.Ring) *MultiArray {
	a.ApplyBinaryOp(NTTMul{}, q, a, baseRing)
	return a
}

// Neg negates the polynomial in-place.
func (a *MultiArray) Neg(baseRing *ring.Ring) *MultiArray {
	// TODO implement.
	return a
}

// -- Vector helpers

// DotProduct performs a coefficient-wise dot product of the vectors and returns the result.
func (v Vector) DotProduct(b Vector, baseRing *ring.Ring) Polynomial {
	out := NewPolynomial(baseRing)
	v.ForEach(func(el Polynomial, i int) {
		el.ApplyBinaryOp(NTTMulAdd{}, b.ElementAtIndex(i), out, baseRing)
	})
	return out
}

// -- Matrix helpers

// MulVecTo performs a matrix vector multiplication and returns the result.
func (m Matrix) MulVec(x Vector, mulAddOp BinaryOperator, baseRing *ring.Ring) Vector {
	out := NewVectorFromSize(m.Rows(), baseRing)
	m.ForEachRow(func(row Vector, i int) {
		out.SetElementAtIndex(i, row.DotProduct(x, baseRing))
	})
	return out
}

// -- Conversion helpers

func (a *MultiArray) AsMatrix() Matrix {
	// assert len(a.Dimensions()) == 2
	return Matrix{a}
}

func (a *MultiArray) AsVector() Vector {
	// assert len(a.Dimensions()) == 1
	return Vector{a}
}
