package math

// Helpers for specialized polynomial constructs

type PolyArray struct {
	*MultiArray
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

// Scale scales the coefficients of the polynomials in-place.
// p_i => c*p_i
func (m PolyArray) Scale(factor uint64) PolyArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).Scale(factor)
	})
	return m
}

// -- Conversion helpers

// AsMatrix converts the representation to a matrix.
func (m *MultiArray) AsMatrix() Matrix {
	// assert len(a.Dimensions()) == 2
	return Matrix{m}
}

// AsVector converts the representation to a vector.
func (m *MultiArray) AsVector() Vector {
	// assert len(a.Dimensions()) == 1
	return Vector{m}
}

// AsPolyArray converts the representation of a multi array of polynomials.
// Allows calling specialized methods.
func (m *MultiArray) AsPolyArray() PolyArray {
	// assert type
	return PolyArray{m}
}
