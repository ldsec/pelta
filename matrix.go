package main

import "github.com/ldsec/lattigo/v2/ring"

// MultiArray represents a multidimensional array of polynomials.
type MultiArray struct {
	dims  []int
	mults []int        // Multiplicants for each coordinate.
	Array []*ring.Poly // Linearized Array.
}

// NewMultiArray constructs a new matrix with the given dimensions.
func NewMultiArray(dims []int, baseRing *ring.Ring) MultiArray {
	totalLength := 1
	for i := 0; i < len(dims); i++ {
		totalLength *= dims[i]
	}
	mults := []int{1}
	acc := 1
	for i := 1; i < len(dims); i++ {
		acc *= dims[i-1]
		mults = append(mults, acc)
	}
	// Initialize the underlying Array to empty polynomials.
	array := make([]*ring.Poly, totalLength)
	for i := 0; i < len(array); i++ {
		array[i] = baseRing.NewPoly()
	}
	return MultiArray{
		dims:  dims,
		mults: mults,
		Array: array,
	}
}

// -- Vector stuff

func NewVectorFromDimensions(dim int, baseRing *ring.Ring) MultiArray {
	return NewMultiArray([]int{dim}, baseRing)
}

func NewVectorFromSlice(elements []*ring.Poly) MultiArray {
	return MultiArray{
		dims:  []int{len(elements)},
		mults: []int{1},
		Array: elements,
	}
}

func (m *MultiArray) IsVector() bool {
	return len(m.dims) == 1
}

// -- Matrix stuff

func NewMatrixFromDimensions(rows int, cols int, baseRing *ring.Ring) MultiArray {
	return NewMultiArray([]int{cols, rows}, baseRing)
}

func NewMatrixFromSlice(array [][]*ring.Poly, baseRing *ring.Ring) MultiArray {
	// TODO: implement!
	return MultiArray{}
}

func (m *MultiArray) IsMatrix() bool {
	return len(m.dims) == 2
}

func (m *MultiArray) MatrixRowSlice(row int) MultiArray {
	indexStart := m.fromCoords([]int{0, row})
	indexEnd := m.fromCoords([]int{0, row + 1})
	rowArray := m.Array[indexStart:indexEnd]
	return MultiArray{Array: rowArray, mults: []int{1}, dims: []int{len(rowArray)}}
}

func (m *MultiArray) MatrixElement(row int, col int) *ring.Poly {
	index := m.fromCoords([]int{col, row})
	return m.Array[index]
}

func (m *MultiArray) MatrixSetRow(row int, v MultiArray) {
	indexStart := m.fromCoords([]int{0, row})
	indexEnd := m.fromCoords([]int{0, row + 1})
	for i := indexStart; i < indexEnd; i++ {
		m.SetElementAtIndex(i, v.Array[i])
	}
}

func (m *MultiArray) MatrixTranspose() MultiArray {
	// TODO: implement!
	return MultiArray{}
}

// --

func (m *MultiArray) Dimensions() []int {
	return m.dims
}

func (m *MultiArray) Length() int {
	return len(m.Array)
}

func (m *MultiArray) ElementAtCoords(coords []int) *ring.Poly {
	index := m.fromCoords(coords)
	return m.Array[index]
}

func (m *MultiArray) ElementAtIndex(index int) *ring.Poly {
	return m.Array[index]
}

func (m *MultiArray) SetElementAtCoords(coords []int, newElement *ring.Poly) {
	index := m.fromCoords(coords)
	m.Array[index] = newElement
}

func (m *MultiArray) SetElementAtIndex(index int, newElement *ring.Poly) {
	m.Array[index] = newElement
}

func (m *MultiArray) Map(f func(*ring.Poly, []int) *ring.Poly) {
	for i := 0; i < len(m.Array); i++ {
		coords := m.toCoords(i)
		m.Array[i] = f(m.Array[i], coords)
	}
}

func (m *MultiArray) ForEach(f func(*ring.Poly, []int)) {
	for i := 0; i < len(m.Array); i++ {
		coords := m.toCoords(i)
		f(m.Array[i], coords)
	}
}

// -- Simple base conversion in-between the coordinate space and the index space.

func (m *MultiArray) fromCoords(coords []int) int {
	acc := 0
	for i := 0; i < len(m.dims); i++ {
		acc += coords[i] * m.mults[i]
	}
	return acc
}

func (m *MultiArray) toCoords(linearIndex int) []int {
	acc := linearIndex
	coords := make([]int, len(m.dims))
	for i := len(m.mults) - 1; i >= 0; i-- {
		qtd := acc / m.mults[i]
		rem := acc % m.mults[i]
		coords[i] = qtd
		acc = rem
	}
	return coords
}
