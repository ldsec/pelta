package math

// -- Polynomial helpers

// -- Multiarray of poly helpers

// NTT converts the vector of polynomials into the NTT space in-place.
// p_i => NTT(p_i)
func (m *MultiArray) NTT() *MultiArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).NTT()
	})
	return m
}

// InvNTT converts the vector of NTT polynomials back into the poly space in-place.
// p_i => InvNTT(p_i)
func (m *MultiArray) InvNTT() *MultiArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).InvNTT()
	})
	return m
}

// Scale scales the coefficients of the polynomials in-place.
// p_i => c*p_i
func (m *MultiArray) Scale(factor uint64) *MultiArray {
	m.ForEach(func(el RingElement, _ []int) {
		el.(Polynomial).Scale(factor)
	})
	return m
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
