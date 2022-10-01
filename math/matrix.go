package math

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Matrix represents a 2-dimensional matrix.
type Matrix struct {
	*MultiArray
}

// NewMatrixFromDimensions constructs an empty matrix with the given dimensions.
func NewMatrixFromDimensions(rows int, cols int, baseRing *ring.Ring) Matrix {
	return Matrix{NewMultiArray([]int{cols, rows}, baseRing)}
}

// NewMatrixFromSlice constructs a matrix from the given slice.
// The slice is assumed to be a list of row vectors.
func NewMatrixFromSlice(array [][]Polynomial, baseRing *ring.Ring) Matrix {
	// Create an empty matrix.
	rows, cols := len(array), len(array[0])
	m := NewMatrixFromDimensions(rows, cols, baseRing)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			m.SetElement(i, j, array[i][j])
		}
	}
	return m
}

// PopulateRows initializes the rows of this multi array using a given function.
func (m Matrix) PopulateRows(f func(int) Vector) Matrix {
	for i := 0; i < m.Rows(); i++ {
		m.SetRow(i, f(i))
	}
	return m
}

// PopulateCols initializes the cols of this multi array using a given function.
func (m Matrix) PopulateCols(f func(int) Vector) Matrix {
	for i := 0; i < m.Cols(); i++ {
		m.SetCol(i, f(i))
	}
	return m
}

// Rows returns the number of rows this matrix has.
func (m Matrix) Rows() int {
	return m.coordMap.dims[1]
}

// Cols returns the number of cols this matrix has.
func (m Matrix) Cols() int {
	return m.coordMap.dims[0]
}

// Dimensions returns the dimensions of this matrix as a row, col pair.
func (m *Matrix) Dimensions() (int, int) {
	return m.Rows(), m.Cols()
}

// Row returns the vector representation of the given row.
func (m Matrix) Row(row int) Vector {
	// assert row < m.Rows()
	indexStart := m.coordMap.FromCoords([]int{0, row})
	indexEnd := m.coordMap.FromCoords([]int{0, row + 1})
	rowArray := m.Array[indexStart:indexEnd]
	return NewVectorFromSlice(rowArray, m.baseRing)
}

// Col returns the vector representation of the given col.
func (m Matrix) Col(col int) Vector {
	// assert col < m.Cols()
	colSlice := make([]Polynomial, m.Rows())
	for i := 0; i < m.Rows(); i++ {
		colSlice[i] = m.Element(i, col)
	}
	return NewVectorFromSlice(colSlice, m.baseRing)
}

// Element returns the element at the given position.
func (m Matrix) Element(row int, col int) Polynomial {
	return m.ElementAtCoords([]int{col, row})
}

// SetElement updates the element at the given position.
func (m Matrix) SetElement(row int, col int, newElement Polynomial) {
	m.SetElementAtCoords([]int{col, row}, newElement)
}

// SetRow updates the row with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m Matrix) SetRow(row int, v Vector) {
	// assert m.Cols() == v.Length()
	m.SetElements([]int{0, row}, []int{0, row + 1}, v.Array)
}

// SetCol updates the col with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m Matrix) SetCol(col int, v Vector) {
	// assert m.Rows() == v.Length()
	for i := 0; i < m.Rows(); i++ {
		m.SetElement(i, col, v.ElementAtIndex(i))
	}
}

// TransposeInPlace transposes the matrix in-place.
func (m Matrix) TransposeInPlace() Matrix {
	// TODO can optimize to use swaps.
	// Shallow copy the elements in.
	old := make([]Polynomial, m.Length())
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

// MapInPlace replaces every cell with the output of the given function in-place.
func (m Matrix) MapInPlace(f func(Polynomial, int, int) Polynomial) Matrix {
	m.MultiArray.MapInPlace(func(el Polynomial, coords []int) Polynomial {
		return f(el, coords[1], coords[0])
	})
	return m
}

// MapRowsInPlace replaces every row with the output of the given function in-place.
func (m Matrix) MapRowsInPlace(f func(Vector, int) Vector) Matrix {
	for i := 0; i < m.Rows(); i++ {
		m.SetRow(i, f(m.Row(i), i))
	}
	return m
}

// ForEach calls the given function with the contents of each cell.
func (m Matrix) ForEach(f func(Polynomial, int, int)) {
	m.MultiArray.ForEach(func(el Polynomial, coords []int) {
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
