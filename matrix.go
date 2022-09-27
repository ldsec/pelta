package main

import "github.com/ldsec/lattigo/v2/ring"

// MultiArray represents a multidimensional array of polynomials.
type MultiArray struct {
	dims  []int
	mults []int        // Multiplicants for each coordinate.
	array []*ring.Poly // Linearized array.
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
	// Initialize the underlying array to empty polynomials.
	array := make([]*ring.Poly, totalLength)
	for i := 0; i < len(array); i++ {
		array[i] = baseRing.NewPoly()
	}
	return MultiArray{
		dims:  dims,
		mults: mults,
		array: array,
	}
}

func (m *MultiArray) IsVector() bool {
	return len(m.dims) == 1
}

func (m *MultiArray) Dimensions() []int {
	return m.dims
}

func (m *MultiArray) Length() int {
	return len(m.array)
}

func (m *MultiArray) ElementAtCoord(coords []int) *ring.Poly {
	index := m.fromCoords(coords)
	return m.array[index]
}

func (m *MultiArray) ElementAtIndex(index int) *ring.Poly {
	return m.array[index]
}

func (m *MultiArray) SetElement(coords []int, newElement *ring.Poly) error {
	index := m.fromCoords(coords)
	m.array[index] = newElement
	return nil
}

func (m *MultiArray) Map(f func(*ring.Poly, []int) *ring.Poly) {
	for i := 0; i < len(m.array); i++ {
		coords := m.toCoords(i)
		m.array[i] = f(m.array[i], coords)
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
