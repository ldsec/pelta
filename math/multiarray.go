package math

import "github.com/ldsec/lattigo/v2/ring"

// MultiArray represents a multidimensional array of polynomials.
type MultiArray struct {
	coordMap CoordMap
	Array    []Polynomial // Linearized array.
}

// NewMultiArray constructs a new empty multi array with the given dimensions.
func NewMultiArray(dims []int, baseRing *ring.Ring) MultiArray {
	totalLength := 1
	for i := 0; i < len(dims); i++ {
		totalLength *= dims[i]
	}
	// Initialize the underlying Array to empty polynomials.
	array := make([]Polynomial, totalLength)
	for i := 0; i < len(array); i++ {
		array[i] = Polynomial{baseRing.NewPoly()}
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
func (m *MultiArray) ElementAtCoords(coords []int) Polynomial {
	index := m.coordMap.FromCoords(coords)
	return m.Array[index]
}

// ElementAtIndex returns the element at the given linear index.
func (m *MultiArray) ElementAtIndex(index int) Polynomial {
	return m.Array[index]
}

// SetElementAtCoords updates the element at the given coordinates.
func (m *MultiArray) SetElementAtCoords(coords []int, newElement Polynomial) {
	index := m.coordMap.FromCoords(coords)
	m.Array[index] = newElement
}

// SetElementAtIndex updates the element at the given linear index.
func (m *MultiArray) SetElementAtIndex(index int, newElement Polynomial) {
	m.Array[index] = newElement
}

// SetElements updates the array within the given coordinates.
func (m *MultiArray) SetElements(startCoords []int, endCoords []int, replacement []Polynomial) {
	startIndex := m.coordMap.FromCoords(startCoords)
	endIndex := m.coordMap.FromCoords(endCoords)
	// assert |endIndex - startIndex| == len(replacement)
	for i := 0; i < m.abs(endIndex-startIndex); i++ {
		m.Array[m.min(startIndex, endIndex)+i] = replacement[i]
	}
}

// MapInPlace maps every cell with as the output given function in-place.
func (m *MultiArray) MapInPlace(f func(Polynomial, []int) Polynomial) {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		// Update the value in-place.
		m.Array[i] = f(m.Array[i], coords)
	}
}

// ForEach calls the given function with the contents of each cell.
func (m *MultiArray) ForEach(f func(Polynomial, []int)) {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		f(m.Array[i], coords)
	}
}

// ApplyUnaryOp applies the given unary operator to every cell.
func (m *MultiArray) ApplyUnaryOp(op UnaryOperator, out *MultiArray, baseRing *ring.Ring) {
	// assert m.Length() == out.Length()
	for i := 0; i < m.Length(); i++ {
		m.Array[i].ApplyUnaryOp(op, out.Array[i], baseRing)
	}
}

// ApplyBinaryOp applies the given binary operator to this matrix with the elements of `r`, element-wise.
func (m *MultiArray) ApplyBinaryOp(op BinaryOperator, r *MultiArray, out *MultiArray, baseRing *ring.Ring) {
	// assert r.Length() == m.Length() == out.Length()
	for i := 0; i < m.Length(); i++ {
		m.Array[i].ApplyBinaryOp(op, r.Array[i], out.Array[i], baseRing)
	}
}

// ApplyReduction reduces the elements into a single one by repeatedly applying the given operator.
// Warning: Does not clear the `out` polynomial.
func (m *MultiArray) ApplyReduction(op BinaryOperator, out Polynomial, baseRing *ring.Ring) {
	// assert m.Length() == out.Length()
	for i := 0; i < m.Length(); i++ {
		out.ApplyBinaryOp(op, m.Array[i], out, baseRing)
	}
}

// DeepCopy returns a deep copy of the multi array, where each element is also copied.
func (m *MultiArray) DeepCopy() MultiArray {
	// Create an exact copy of the contents.
	array := make([]Polynomial, m.Length())
	for i := 0; i < len(array); i++ {
		array[i] = m.Array[i].DeepCopy()
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
