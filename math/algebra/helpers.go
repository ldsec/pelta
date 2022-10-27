package algebra

import "fmt"

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
