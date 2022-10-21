package math

import "fmt"

// Matrix represents a 2-dimensional matrix.
type Matrix struct {
	*MultiArray
}

// NewMatrixFromDimensions constructs an empty matrix with the given dimensions.
func NewMatrixFromDimensions(rows int, cols int) Matrix {
	return Matrix{NewMultiArray([]int{cols, rows})}
}

// NewMatrixFromSlice constructs a matrix from the given slice.
// The slice is assumed to be a list of row vectors.
func NewMatrixFromSlice(array [][]RingElement) Matrix {
	// Create an empty matrix.
	rows, cols := len(array), len(array[0])
	m := NewMatrixFromDimensions(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			m.SetElement(i, j, array[i][j])
		}
	}
	return m
}

// Populate initializes the elements of this multi array using a given function.
func (m Matrix) Populate(f func(int, int) RingElement) Matrix {
	m.MultiArray.Populate(func(coords []int) RingElement {
		col, row := coords[0], coords[1]
		return f(row, col)
	})
	return m
}

// PopulateRows initializes the rows of this multi array using a given function.
func (m Matrix) PopulateRows(f func(int) Vector) Matrix {
	for i := 0; i < m.Rows(); i++ {
		newRow := f(i)
		if newRow.Length() != m.Cols() {
			panic(fmt.Sprintf("Matrix.PopulateRows: Returned invalid row %d != %d", newRow.Length(), m.Cols()))
		}
		m.SetRow(i, newRow)
	}
	return m
}

// PopulateCols initializes the cols of this multi array using a given function.
func (m Matrix) PopulateCols(f func(int) Vector) Matrix {
	for i := 0; i < m.Cols(); i++ {
		newCol := f(i)
		if newCol.Length() != m.Rows() {
			panic(fmt.Sprintf("Matrix.PopulateCols: Returned invalid col %d != %d", newCol.Length(), m.Rows()))
		}
		m.SetCol(i, newCol)
	}
	return m
}

// Rows returns the number of rows this matrix has.
func (m Matrix) Rows() int {
	return m.coordMap.Dims[1]
}

// Cols returns the number of cols this matrix has.
func (m Matrix) Cols() int {
	return m.coordMap.Dims[0]
}

// Dimensions returns the dimensions of this matrix as a row, col pair.
func (m Matrix) Dimensions() (int, int) {
	return m.Rows(), m.Cols()
}

// Row returns the vector representation of the given row.
func (m Matrix) Row(row int) Vector {
	// assert row < m.Rows()
	if row >= m.Rows() {
		panic(fmt.Sprintf("Matrix.Row: Invalid row %d >= %d", row, m.Rows()))
	}
	indexStart := m.coordMap.FromCoords([]int{0, row})
	indexEnd := m.coordMap.FromCoords([]int{0, row + 1})
	rowArray := m.Array[indexStart:indexEnd]
	return NewVectorFromSlice(rowArray)
}

// Col returns the vector representation of the given col.
func (m Matrix) Col(col int) Vector {
	// assert col < m.Cols()
	if col >= m.Cols() {
		panic(fmt.Sprintf("Matrix.Col: Invalid col %d >= %d", col, m.Cols()))
	}
	colSlice := make([]RingElement, m.Rows())
	for i := 0; i < m.Rows(); i++ {
		colSlice[i] = m.Element(i, col)
	}
	return NewVectorFromSlice(colSlice)
}

// Element returns the element at the given position.
func (m Matrix) Element(row int, col int) RingElement {
	return m.ElementAtCoords([]int{col, row})
}

// SetElement updates the element at the given position.
func (m Matrix) SetElement(row int, col int, newElement RingElement) {
	m.SetElementAtCoords([]int{col, row}, newElement)
}

// SetRow updates the row with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m Matrix) SetRow(row int, v Vector) {
	// assert m.Cols() == v.Length()
	if m.Cols() != v.Length() {
		panic(fmt.Sprintf("Matrix.SetRow: Lengths do not match %d != %d", m.Cols(), v.Length()))
	}
	m.SetElements([]int{0, row}, []int{0, row + 1}, v.Array)
}

// SetCol updates the col with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m Matrix) SetCol(col int, v Vector) {
	// assert m.Rows() == v.Length()
	if m.Rows() != v.Length() {
		panic(fmt.Sprintf("Matrix.SetCol: Lengths do not match %d != %d", m.Rows(), v.Length()))
	}
	for i := 0; i < m.Rows(); i++ {
		m.SetElement(i, col, v.ElementAtIndex(i))
	}
}

// Transpose transposes the matrix in-place.
func (m Matrix) Transpose() Matrix {
	// TODO can optimize to use swaps.
	// Shallow copy the elements in.
	old := make([]RingElement, m.Length())
	for i := 0; i < m.Length(); i++ {
		old[i] = m.Array[i]
	}
	// Swap the dimensions, keep the old coordinate map.
	oldCordMaap := m.coordMap.Copied()
	m.coordMap = m.coordMap.Reversed()
	// Copy the array in, with A^T_ij = A_ji
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			// Note that row is always the second dimension.
			oldIndex := oldCordMaap.FromCoords([]int{i, j})
			m.SetElement(i, j, old[oldIndex])
		}
	}
	return m
}

// Map replaces every cell with the output of the given function in-place.
func (m Matrix) Map(f func(RingElement, int, int) RingElement) Matrix {
	return m.MultiArray.Map(func(el RingElement, coords []int) RingElement {
		return f(el, coords[1], coords[0])
	}).AsMatrix()
}

// MapRows replaces every row with the output of the given function in-place.
func (m Matrix) MapRows(f func(Vector, int) Vector) Matrix {
	for i := 0; i < m.Rows(); i++ {
		newRow := f(m.Row(i), i)
		if newRow.Length() >= m.Cols() {
			panic(fmt.Sprintf("Matrix.MapRows: Returned invalid row %d >= %d", newRow.Length(), m.Cols()))
		}
		m.SetRow(i, newRow)
	}
	return m
}

// ForEach calls the given function with the contents of each cell.
func (m Matrix) ForEach(f func(RingElement, int, int)) {
	m.MultiArray.ForEach(func(el RingElement, coords []int) {
		f(el, coords[1], coords[0])
	})
}

// ForEachRow calls the given function for each row vector.
func (m Matrix) ForEachRow(f func(Vector, int)) {
	for i := 0; i < m.Rows(); i++ {
		targetRow := m.Row(i)
		f(targetRow, i)
	}
}

// ForEachCol calls the given function for each col vector.
func (m Matrix) ForEachCol(f func(Vector, int)) {
	for i := 0; i < m.Cols(); i++ {
		targetCol := m.Col(i)
		f(targetCol, i)
	}
}

// MulVec performs a matrix vector multiplication and returns the result.
func (m Matrix) MulVec(x Vector) Vector {
	if m.Cols() != x.Length() {
		panic(fmt.Sprintf("Matrix.MulVec: Cannot multiply (%d, %d) with %d", m.Rows(), m.Cols(), x.Length()))
	}
	out := NewVectorFromSize(m.Rows())
	m.ForEachRow(func(row Vector, i int) {
		out.SetElementAtIndex(i, row.Dot(x))
	})
	return out
}

// All returns true if all the elements return true on the given predicate.
func (m Matrix) All(pred func(RingElement, int, int) bool) bool {
	return m.MultiArray.All(func(el RingElement, coords []int) bool {
		return pred(el, coords[1], coords[0])
	})
}

// AllRows returns true if all the rows return true on the given predicate.
func (m Matrix) AllRows(pred func(Vector, int) bool) bool {
	for i := 0; i < m.Rows(); i++ {
		if !pred(m.Row(i), i) {
			return false
		}
	}
	return true
}

// AllCols returns true if all the cols return true on the given predicate.
func (m Matrix) AllCols(pred func(Vector, int) bool) bool {
	for i := 0; i < m.Cols(); i++ {
		if !pred(m.Col(i), i) {
			return false
		}
	}
	return true
}

// Any returns true if some element returns true on the given predicate.
func (m Matrix) Any(pred func(RingElement, int, int) bool) bool {
	return m.MultiArray.Any(func(el RingElement, coords []int) bool {
		return pred(el, coords[1], coords[0])
	})
}

// AnyRows returns true if some row returns true on the given predicate.
func (m Matrix) AnyRows(pred func(Vector, int) bool) bool {
	for i := 0; i < m.Rows(); i++ {
		if pred(m.Row(i), i) {
			return true
		}
	}
	return false
}

// AnyCols returns true if some col returns true on the given predicate.
func (m Matrix) AnyCols(pred func(Vector, int) bool) bool {
	for i := 0; i < m.Cols(); i++ {
		if pred(m.Col(i), i) {
			return true
		}
	}
	return false
}
