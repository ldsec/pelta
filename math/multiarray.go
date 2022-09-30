package math

import "github.com/ldsec/lattigo/v2/ring"

// MultiArray represents a multidimensional array of polynomials.
type MultiArray struct {
	coordMap CoordMap
	Array    []*ring.Poly // Linearized array.
}

// NewMultiArray constructs a new empty multi array with the given dimensions.
func NewMultiArray(dims []int, baseRing *ring.Ring) MultiArray {
	totalLength := 1
	for i := 0; i < len(dims); i++ {
		totalLength *= dims[i]
	}
	// Initialize the underlying Array to empty polynomials.
	array := make([]*ring.Poly, totalLength)
	for i := 0; i < len(array); i++ {
		array[i] = baseRing.NewPoly()
	}
	return MultiArray{
		coordMap: NewCoordMap(dims),
		Array:    array,
	}
}

// Dimensions returns the dimensions of this multi array.
func (m *MultiArray) Dimensions() []int {
	return m.coordMap.dims
}

// Length returns the total number of cells in this multi array.
func (m *MultiArray) Length() int {
	return len(m.Array)
}

// ElementAtCoords returns the element at the given coordinates.
func (m *MultiArray) ElementAtCoords(coords []int) *ring.Poly {
	index := m.coordMap.FromCoords(coords)
	return m.Array[index]
}

// ElementAtIndex returns the element at the given linear index.
func (m *MultiArray) ElementAtIndex(index int) *ring.Poly {
	return m.Array[index]
}

// SetElementAtCoords updates the element at the given coordinates.
func (m *MultiArray) SetElementAtCoords(coords []int, newElement *ring.Poly) {
	index := m.coordMap.FromCoords(coords)
	m.Array[index] = newElement
}

// SetElementAtIndex updates the element at the given linear index.
func (m *MultiArray) SetElementAtIndex(index int, newElement *ring.Poly) {
	m.Array[index] = newElement
}

// SetElements updates the array within the given coordinates.
func (m *MultiArray) SetElements(startCoords []int, endCoords []int, replacement []*ring.Poly) {
	startIndex := m.coordMap.FromCoords(startCoords)
	endIndex := m.coordMap.FromCoords(endCoords)
	// assert |endIndex - startIndex| == len(replacement)
	for i := 0; i < m.abs(endIndex-startIndex); i++ {
		m.Array[m.min(startIndex, endIndex)+i] = replacement[i]
	}
}

// MapInPlace maps every cell with as the output given function in-place.
func (m *MultiArray) MapInPlace(f func(*ring.Poly, []int) *ring.Poly) {
	for i := 0; i < len(m.Array); i++ {
		coords := m.coordMap.ToCoords(i)
		// Update the value in-place.
		m.Array[i] = f(m.Array[i], coords)
	}
}

// ForEach calls the given function with the contents of each cell.
func (m *MultiArray) ForEach(f func(*ring.Poly, []int)) {
	for i := 0; i < len(m.Array); i++ {
		coords := m.coordMap.ToCoords(i)
		f(m.Array[i], coords)
	}
}

// ApplyUnaryOp applies the given unary operator to every cell in-place.
func (m *MultiArray) ApplyUnaryOp(op UnaryOperator, baseRing *ring.Ring) {
	for i := 0; i < len(m.Array); i++ {
		// Update the value in-place.
		op.Apply(m.Array[i], m.Array[i], baseRing)
	}
}

// DeepCopy returns a deep copy of the multi array, where each element is also copied.
func (m *MultiArray) DeepCopy() MultiArray {
	// Create an exact copy of the contents.
	array := make([]*ring.Poly, m.Length())
	for i := 0; i < len(array); i++ {
		array[i] = m.Array[i].CopyNew()
	}
	return MultiArray{
		coordMap: m.coordMap.Copied(),
		Array:    array,
	}
}

// ShallowCopy returns a shallow copy of the multi array.
func (m *MultiArray) ShallowCopy() MultiArray {
	return MultiArray{
		coordMap: m.coordMap.Copied(),
		Array:    m.Array,
	}
}
