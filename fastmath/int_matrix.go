package fastmath

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type IntMatrix struct {
	numRows  int
	numCols  int
	rows     []*IntVec
	baseRing *ring.Ring

	CachedTranspose *IntMatrix
}

func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) *IntMatrix {
	rows := make([]*IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = NewIntVec(numCols, baseRing)
	}
	return &IntMatrix{numRows, numCols, rows, baseRing, nil}
}

// NewIdIntMatrix returns an n by n identity matrix.
func NewIdIntMatrix(numRows int, baseRing *ring.Ring) *IntMatrix {
	m := NewIntMatrix(numRows, numRows, baseRing)
	for i, r := range m.RowsView() {
		r.SetForce(i, 1)
	}
	return m
}

func NewIntMatrixFromRows(rows []*IntVec, baseRing *ring.Ring) *IntMatrix {
	return &IntMatrix{len(rows), rows[0].Size(), rows, baseRing, nil}
}

func NewIntMatrixFromSlice(elems [][]uint64, baseRing *ring.Ring) *IntMatrix {
	numRows := len(elems)
	numCols := len(elems[0])
	m := NewIntMatrix(numRows, numCols, baseRing)
	m.PopulateRows(func(i int) *IntVec {
		return NewIntVecFromSlice(elems[i], baseRing)
	})
	return m
}

func NewIntMatrixFromCoeffSlice(elems [][]Coeff, baseRing *ring.Ring) *IntMatrix {
	numRows := len(elems)
	numCols := len(elems[0])
	m := NewIntMatrix(numRows, numCols, baseRing)
	m.PopulateRows(func(i int) *IntVec {
		return NewIntVecFromCoeffSlice(elems[i], baseRing)
	})
	return m
}

func (m *IntMatrix) BaseRing() *ring.Ring {
	return m.baseRing
}

// Rows returns the number of rows.
func (m *IntMatrix) Rows() int {
	return m.numRows
}

// Cols returns the number of cols.
func (m *IntMatrix) Cols() int {
	return m.numCols
}

// RowsView returns a list of references to the rows of this matrix.
func (m *IntMatrix) RowsView() []*IntVec {
	return m.rows
}

// RowView returns a reference to the i-th row.
func (m *IntMatrix) RowView(i int) *IntVec {
	return m.rows[i]
}

// SubsectionCopy returns a subsection (copied) of this matrix.
func (m *IntMatrix) SubsectionCopy(rowStart, rowEnd int, colStart, colEnd int) *IntMatrix {
	subMatrix := NewIntMatrix(rowEnd-rowStart, colEnd-colStart, m.baseRing)
	for i := 0; i < rowEnd-rowStart; i++ {
		subMatrix.SetRow(i, m.RowView(rowStart+i).SliceCopy(colStart, colEnd))
	}
	return subMatrix
}

// ColCopy returns a copy of the i-th col.
func (m *IntMatrix) ColCopy(i int) *IntVec {
	colVec := NewIntVec(m.Rows(), m.baseRing)
	for j, row := range m.rows {
		colVec.SetCoeff(j, row.GetCoeff(i))
	}
	return colVec
}

// GetCoeff returns the element at the given coordinates.
func (m *IntMatrix) GetCoeff(row, col int) Coeff {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Get indices incorrect")
	}
	return m.rows[row].GetCoeff(col)
}

// Get returns the element at the given coordinates.
func (m *IntMatrix) GetLevel(row, col, level int) uint64 {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Get indices incorrect")
	}
	return m.rows[row].GetLevel(col, level)
}

// SetForce updates the given element of this matrix.
func (m *IntMatrix) SetForce(row, col int, newValue uint64) {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	m.rows[row].SetForce(col, newValue)
}

// SetCoeff updates the given element of this matrix.
func (m *IntMatrix) SetCoeff(row, col int, newCoeff Coeff) {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	m.rows[row].SetCoeff(col, newCoeff)
}

// SetRow updates the given row of this matrix.
func (m *IntMatrix) SetRow(row int, newRow *IntVec) {
	if row >= m.Rows() || newRow.Size() != m.Cols() {
		panic("IntMatrix.SetRows index or size incorrect")
	}
	m.rows[row] = newRow
}

// SetCol updates the given col of this matrix.
func (m *IntMatrix) SetCol(col int, newCol *IntVec) {
	if col >= m.Cols() || newCol.Size() != m.Rows() {
		panic("IntMatrix.SetCols index or size incorrect")
	}
	for i, row := range m.rows {
		row.SetCoeff(col, newCol.GetCoeff(i))
	}
}

// AppendRow appends a row into the vector.
func (m *IntMatrix) AppendRow(v *IntVec) {
	if m.Cols() != v.Size() {
		panic("IntMatrix.AppendRow cannot append row, invalid size")
	}
	m.rows = append(m.rows, v)
	m.numRows += 1
}

// ExtendRows concatenates the matrices vertically.
func (m *IntMatrix) ExtendRows(b *IntMatrix) *IntMatrix {
	if m.Cols() != b.Cols() {
		panic("IntMatrix.ExtendRows cannot extend, invalid col size")
	}
	m.rows = append(m.rows, b.rows...)
	m.numRows += b.Rows()
	return m
}

// ExtendCols concatenates the matrices horizontally.
func (m *IntMatrix) ExtendCols(b *IntMatrix) *IntMatrix {
	if m.Rows() != b.Rows() {
		panic("IntMatrix.ExtendCols cannot extend, invalid row size")
	}
	for i := 0; i < m.Rows(); i++ {
		m.rows[i].Append(b.RowView(i))
	}
	m.numCols += b.Cols()
	return m
}

// Populate is used to initialize the elements of this matrix.
func (m *IntMatrix) Populate(f func(int, int) uint64) {
	for i, row := range m.rows {
		row.Populate(
			func(j int) uint64 {
				return f(i, j)
			})
	}
}

// PopulateRows is used to initialize the rows of this matrix.
func (m *IntMatrix) PopulateRows(f func(int) *IntVec) {
	for i := 0; i < m.Rows(); i++ {
		m.rows[i] = f(i)
	}
}

// Copy returns a copy of this matrix.
func (m *IntMatrix) Copy() *IntMatrix {
	rows := make([]*IntVec, m.numRows)
	for i, row := range m.rows {
		rows[i] = row.Copy()
	}
	return &IntMatrix{m.numRows, m.numCols, rows, m.baseRing, nil}
}

// String returns a string representation of the matrix.
func (m *IntMatrix) String() string {
	s := fmt.Sprintf("IntMatrix[%d,%d]{\n", m.Rows(), m.Cols())
	for _, row := range m.rows {
		s += "\t" + row.String() + "\n"
	}
	return s + ", ...}"
}

// SizeString returns a string representation of the matrix's dimensions.
func (m *IntMatrix) SizeString() string {
	s := fmt.Sprintf("IntMatrix[%d,%d]", m.Rows(), m.Cols())
	return s
}

// RebaseRowsLossless rebases each row.
func (m *IntMatrix) RebaseRowsLossless(newRing RingParams) *IntMatrix {
	return logging.LogShortExecution("IntMatrix.RebaseRowsLossless", "rebasing", func() interface{} {
		m.baseRing = newRing.BaseRing
		for _, row := range m.rows {
			row.RebaseLossless(newRing)
		}
		return m
	}).(*IntMatrix)
}

// Scale scales the matrix by the given amount.
func (m *IntMatrix) Scale(factor uint64) *IntMatrix {
	for _, row := range m.rows {
		row.Scale(factor)
	}
	return m
}

func (m *IntMatrix) ScaleCoeff(factors Coeff) *IntMatrix {
	for _, row := range m.rows {
		row.ScaleCoeff(factors)
	}
	return m
}

// Neg negates this matrix.
func (m *IntMatrix) Neg() *IntMatrix {
	for _, row := range m.rows {
		row.Neg()
	}
	return m
}

// Transposed returns the transposed version of this matrix.
func (m *IntMatrix) Transposed() *IntMatrix {
	return logging.LogShortExecution("IntMatrix.Transposed", m.SizeString(), func() interface{} {
		At := make([]Coeff, m.Cols()*m.Rows())
		for row := 0; row < m.Rows(); row++ {
			for col := 0; col < m.Cols(); col++ {
				index := col*m.Rows() + row
				At[index] = m.GetCoeff(row, col)
			}
		}
		mt := NewIntMatrix(m.Cols(), m.Rows(), m.BaseRing())
		for row := 0; row < mt.Rows(); row++ {
			for col := 0; col < mt.Cols(); col++ {
				mt.SetCoeff(row, col, At[row*mt.Cols()+col])
			}
		}
		return mt
	}).(*IntMatrix)
}

// Hadamard performs coefficient-wise multiplication.
func (m *IntMatrix) Hadamard(b *IntMatrix) *IntMatrix {
	for i, r := range m.rows {
		r.Hadamard(b.RowView(i))
	}
	return m
}

// Max returns the largest element of the matrix.
func (m *IntMatrix) Max() uint64 {
	max := m.rows[0].Max()
	for _, r := range m.rows {
		c := r.Max()
		if c > max {
			max = c
		}
	}
	return max
}

// Eq returns true iff two matrices are equal in their elements.
func (m *IntMatrix) Eq(b *IntMatrix) bool {
	if m.Rows() != b.Rows() || m.Cols() != b.Cols() {
		return false
	}
	for i, row := range m.rows {
		bRow := b.RowView(i)
		if !row.Eq(bRow) {
			return false
		}
	}
	return true
}

// MulVec performs a matrix-vector multiplication.
func (m *IntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	return logging.LogShortExecution("IntMatrix.MulVec", fmt.Sprintf("%s * [%d]", m.SizeString(), v.Size()), func() interface{} {
		out := NewIntVec(m.Rows(), m.baseRing)
		for i, row := range m.rows {
			dotResult := row.Dot(v)
			out.SetCoeff(i, dotResult)
		}
		return out
	}).(*IntVec)
}

// MulMat performs a matrix-matrix multiplication.
func (m *IntMatrix) MulMat(b *IntMatrix) *IntMatrix {
	if m.Cols() != b.Rows() {
		panic("IntMatrix.MulMat sizes incorrect")
	}
	out := NewIntMatrix(m.Rows(), b.Cols(), m.baseRing)
	for i := 0; i < out.Rows(); i++ {
		for j := 0; j < out.Cols(); j++ {
			p := m.RowView(i)
			q := b.ColCopy(j)
			dotResult := q.Dot(p)
			out.SetCoeff(i, j, dotResult)
		}
	}
	return out
}
