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
func (m *MultiArray) Add(q *MultiArray) *MultiArray {
	m.ApplyBinaryOp(SimpleAdd{}, q, m)
	return m
}

// Sum sums the polynomials in the array, returning the result.
func (m *MultiArray) Sum() Polynomial {
	out := NewPolynomial(m.baseRing)
	m.ApplyReduction(SimpleAdd{}, out)
	return out
}

// Product takes the coefficient-wise product of all the elements in the array, returning the result.
func (m *MultiArray) Product() Polynomial {
	out := NewPolynomial(m.baseRing)
	m.ApplyReduction(NTTMul{}, out)
	return out
}

// NTT converts the vector into the NTT space in-place.
// p => p := NTT(p)
func (m *MultiArray) NTT() *MultiArray {
	m.ApplyUnaryOp(NTT{}, m)
	return m
}

// InvNTT converts the elements back into the poly space in-place.
// p => p := InvNTT(p)
func (m *MultiArray) InvNTT() *MultiArray {
	m.ApplyUnaryOp(InvNTT{}, m)
	return m
}

// Scale scales the coefficients of the polynomials in-place.
// p => p := c*p
func (m *MultiArray) Scale(factor uint64) *MultiArray {
	m.ApplyUnaryOp(Scale{factor}, m)
	return m
}

// Mul multiplies the elements coefficient-wise in-place.
// p, q => p *= q
func (m *MultiArray) Mul(q *MultiArray) *MultiArray {
	m.ApplyBinaryOp(NTTMul{}, q, m)
	return m
}

// Neg negates the polynomial in-place.
func (m *MultiArray) Neg() *MultiArray {
	// TODO implement.
	return m
}

// -- Vector helpers

// DotProduct performs a coefficient-wise dot product of the vectors and returns the result.
func (v Vector) DotProduct(b Vector) Polynomial {
	out := NewPolynomial(v.baseRing)
	v.ForEach(func(el Polynomial, i int) {
		el.ApplyBinaryOp(NTTMulAdd{}, b.ElementAtIndex(i), out, v.baseRing)
	})
	return out
}

// -- Matrix helpers

// MulVecTo performs a matrix vector multiplication and returns the result.
func (m Matrix) MulVec(x Vector) Vector {
	out := NewVectorFromSize(m.Rows(), m.baseRing)
	m.ForEachRow(func(row Vector, i int) {
		out.SetElementAtIndex(i, row.DotProduct(x))
	})
	return out
}

// -- Conversion helpers

func (m *MultiArray) AsMatrix() Matrix {
	// assert len(a.Dimensions()) == 2
	return Matrix{m}
}

func (m *MultiArray) AsVector() Vector {
	// assert len(a.Dimensions()) == 1
	return Vector{m}
}
