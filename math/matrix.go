package math

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Matrix represents a 2-dimensional matrix.
type Matrix struct {
	MultiArray
}

// NewMatrixFromDimensions constructs an empty matrix with the given dimensions.
func NewMatrixFromDimensions(rows int, cols int, baseRing *ring.Ring) Matrix {
	return Matrix{NewMultiArray([]int{cols, rows}, baseRing)}
}

// NewMatrixFromSlice constructs a matrix from the given slice.
// The slice is assumed to be a list of row vectors.
func NewMatrixFromSlice(array [][]*ring.Poly, baseRing *ring.Ring) Matrix {
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

// Rows returns the number of rows this matrix has.
func (m *Matrix) Rows() int {
	return m.coordMap.dims[1]
}

// Cols returns the number of cols this matrix has.
func (m *Matrix) Cols() int {
	return m.coordMap.dims[0]
}

// Dimensions returns the dimensions of this matrix as a row, col pair.
func (m *Matrix) Dimensions() (int, int) {
	return m.Rows(), m.Cols()
}

// Row returns the vector representation of the given row.
func (m *Matrix) Row(row int) Vector {
	// assert row < m.Rows()
	indexStart := m.coordMap.FromCoords([]int{0, row})
	indexEnd := m.coordMap.FromCoords([]int{0, row + 1})
	rowArray := m.Array[indexStart:indexEnd]
	return NewVectorFromSlice(rowArray)
}

// Col returns the vector representation of the given col.
func (m *Matrix) Col(col int) Vector {
	// assert col < m.Cols()
	colSlice := make([]*ring.Poly, m.Rows())
	for i := 0; i < m.Rows(); i++ {
		colSlice[i] = m.Element(i, col)
	}
	return NewVectorFromSlice(colSlice)
}

// Element returns the element at the given position.
func (m *Matrix) Element(row int, col int) *ring.Poly {
	return m.ElementAtCoords([]int{col, row})
}

// SetElement updates the element at the given position.
func (m *Matrix) SetElement(row int, col int, newElement *ring.Poly) {
	m.SetElementAtCoords([]int{col, row}, newElement)
}

// SetRow updates the row with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m *Matrix) SetRow(row int, v *Vector) {
	// assert m.Cols() == v.Length()
	m.SetElements([]int{0, row}, []int{0, row + 1}, v.Array)
}

// SetCol updates the col with the given vector.
// Warning: Does not copy the underlying array of the elements!
func (m *Matrix) SetCol(col int, v *Vector) {
	// assert m.Rows() == v.Length()
	for i := 0; i < m.Rows(); i++ {
		m.SetElement(i, col, v.Array[i])
	}
}

// TransposeInPlace transposes the matrix in-place.
func (m *Matrix) TransposeInPlace() {
	// TODO can optimize to use swaps.
	// Shallow copy the elements in.
	old := make([]*ring.Poly, m.Length())
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
}

// MapInPlace replaces every cell with the output of the given function in-place.
func (m *Matrix) MapInPlace(f func(*ring.Poly, int, int) *ring.Poly) {
	m.MultiArray.MapInPlace(func(el *ring.Poly, coords []int) *ring.Poly {
		return f(el, coords[1], coords[0])
	})
}

// ForEach calls the given function with the contents of each cell.
func (m *Matrix) ForEach(f func(*ring.Poly, int, int)) {
	m.MultiArray.ForEach(func(el *ring.Poly, coords []int) {
		f(el, coords[1], coords[0])
	})
}

// MulVecTo performs a matrix vector multiplication and adds the result to the given `out` vector.
func (m *Matrix) MulVecTo(x *Vector, mulAddOp BinaryOperator, out *Vector, baseRing *ring.Ring) {
	outputDim := m.Rows()
	for i := 0; i < outputDim; i++ {
		targetRow := m.Row(i)
		targetRow.DotProductTo(x, mulAddOp, out.ElementAtIndex(i), baseRing)
	}
}

// MulVec performs a matrix vector multiplication and returns the result.
func (m *Matrix) MulVec(x *Vector, mulAddOp BinaryOperator, baseRing *ring.Ring) Vector {
	outputDim := m.Rows()
	y := NewVectorFromDimensions(outputDim, baseRing)
	m.MulVecTo(x, mulAddOp, &y, baseRing)
	return y
}
