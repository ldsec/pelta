package algebra

import (
	"fmt"
)

// MultiArray represents a multidimensional array of elements.
type MultiArray struct {
	coordMap CoordMap
	Array    []Element // Linearized array.
}

// NewMultiArray constructs a new empty multi array with the given dimensions.
func NewMultiArray(dims []int) *MultiArray {
	totalLength := 1
	for i := 0; i < len(dims); i++ {
		totalLength *= dims[i]
	}
	// Initialize the underlying Array to nil values.
	array := make([]Element, totalLength)
	for i := 0; i < len(array); i++ {
		array[i] = nil
	}
	a := MultiArray{
		coordMap: NewCoordMap(dims),
		Array:    array,
	}
	return &a
}

// Populate initializes the cells of this multi array using a given function.
func (m *MultiArray) Populate(f func([]int) Element) *MultiArray {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		m.Array[i] = f(coords)
	}
	return m
}

// Dimensions returns the dimensions of this multi array.
func (m *MultiArray) Dimensions() []int {
	return m.coordMap.Dims
}

// Length returns the total number of cells in this multi array.
func (m *MultiArray) Length() int {
	return len(m.Array)
}

// ElementAtCoords returns the element at the given coordinates.
func (m *MultiArray) ElementAtCoords(coords []int) Element {
	index := m.coordMap.FromCoords(coords)
	return m.ElementAtIndex(index)
}

// ElementAtIndex returns the element at the given linear index.
func (m *MultiArray) ElementAtIndex(index int) Element {
	if index >= m.Length() {
		panic(fmt.Sprintf("MultiArray.ElementAtIndex: Index out of bounds %d >= %d", index, m.Length()))
	}
	return m.Array[index]
}

// SetElementAtCoords updates the element at the given coordinates.
func (m *MultiArray) SetElementAtCoords(coords []int, newElement Element) {
	index := m.coordMap.FromCoords(coords)
	m.SetElementAtIndex(index, newElement)
}

// SetElementAtIndex updates the element at the given linear index.
func (m *MultiArray) SetElementAtIndex(index int, newElement Element) {
	if index >= m.Length() {
		panic(fmt.Sprintf("MultiArray.SetElementAtIndex: Index out of bounds %d >= %d", index, m.Length()))
	}
	m.Array[index] = newElement
}

// SetElements updates the array within the given coordinates. End exclusive.
func (m *MultiArray) SetElements(startCoords []int, endCoords []int, replacement []Element) {
	// assert startCoords < endCoords
	startIndex := m.coordMap.FromCoords(startCoords)
	endIndex := m.coordMap.FromCoords(endCoords)
	// assert abs(endIndex - startIndex) == len(replacement)
	if abs(endIndex-startIndex) != len(replacement) {
		panic(fmt.Sprintf("MultiArray.SetElements: Replacement does not fit exactly %d != %d",
			len(replacement), abs(endIndex-startIndex)))
	}
	for i := 0; i < abs(endIndex-startIndex); i++ {
		m.Array[minInt(startIndex, endIndex)+i] = replacement[i]
	}
}

// Map maps every cell with as the output given function in-place.
func (m *MultiArray) Map(f func(Element, []int) Element) *MultiArray {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		// Update the value in-place.
		m.Array[i] = f(m.Array[i], coords)
	}
	return m
}

// ForEach calls the given function with the contents of each cell.
func (m *MultiArray) ForEach(f func(Element, []int)) {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		f(m.Array[i], coords)
	}
}

// Copy returns a deep copy of the multi array, where each element is also copied.
func (m *MultiArray) Copy() *MultiArray {
	// Create an exact copy of the contents.
	array := make([]Element, m.Length())
	for i := 0; i < len(array); i++ {
		array[i] = m.Array[i].Copy()
	}
	return &MultiArray{
		coordMap: m.coordMap.Copied(),
		Array:    array,
	}
}

// Scale scales every element by the given factor.
func (m *MultiArray) Scale(factor uint64) *MultiArray {
	for i := 0; i < m.Length(); i++ {
		m.ElementAtIndex(i).Scale(factor)
	}
	return m
}

// Add adds two arrays coefficient-wise in-place.
// p, q => p += q
func (m *MultiArray) Add(q *MultiArray) *MultiArray {
	if m.Length() != q.Length() {
		panic(fmt.Sprintf("MultiArray.Add: Different sizes %d != %d", m.Length(), q.Length()))
	}
	for i := 0; i < m.Length(); i++ {
		m.ElementAtIndex(i).Add(q.ElementAtIndex(i))
	}
	return m
}

// Hadamard multiplies the elements coefficient-wise in-place.
func (m *MultiArray) Hadamard(q *MultiArray) *MultiArray {
	if m.Length() != q.Length() {
		panic(fmt.Sprintf("MultiArray.Hadamard: Different sizes %d != %d", m.Length(), q.Length()))
	}
	for i := 0; i < m.Length(); i++ {
		m.ElementAtIndex(i).Mul(q.ElementAtIndex(i))
	}
	return m
}

// Mul multiplies the elements with the given element in-place.
// p, q => p *= q
func (m *MultiArray) Mul(q Element) *MultiArray {
	m.ForEach(func(el Element, _ []int) {
		el.Mul(q)
	})
	return m
}

// Neg negates the elements in-place.
// p_i => -p_i
func (m *MultiArray) Neg() *MultiArray {
	m.ForEach(func(el Element, _ []int) {
		el.Neg()
	})
	return m
}

// Sum sums the polynomials in the array, returning the result.
func (m *MultiArray) Sum() Element {
	out := m.Array[0].Copy().Zero()
	m.ForEach(func(el Element, _ []int) {
		out.Add(el)
	})
	return out
}

// Product takes the coefficient wise product of all the elements in the array, returning the result.
func (m *MultiArray) Product() Element {
	out := m.Array[0].Copy().One()
	m.ForEach(func(el Element, _ []int) {
		out.Mul(el)
	})
	return out
}

// All returns true if all the elements return true on the given predicate.
func (m *MultiArray) All(pred func(Element, []int) bool) bool {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		if !pred(m.ElementAtIndex(i), coords) {
			return false
		}
	}
	return true
}

// Any returns true if some element returns true on the given predicate.
func (m *MultiArray) Any(pred func(Element, []int) bool) bool {
	for i := 0; i < m.Length(); i++ {
		coords := m.coordMap.ToCoords(i)
		if pred(m.ElementAtIndex(i), coords) {
			return true
		}
	}
	return false
}

// Eq returns true if all the elements are the same in the same positions.
func (m *MultiArray) Eq(other *MultiArray) bool {
	return m.All(func(el Element, coords []int) bool {
		return other.ElementAtCoords(coords).Eq(el)
	})
}