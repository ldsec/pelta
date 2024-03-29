package fastmath

import (
	"fmt"
	"math/big"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// IntMatrix represents a matrix of integers.
type IntMatrix struct {
	numRows      int
	numCols      int
	rows         []*IntVec
	unrebasedRef *IntMatrix
	rebaseRing   *RingParams
	baseRing     *ring.Ring
}

// NewIntMatrix returns an empty int matrix of given size.
func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) *IntMatrix {
	rows := make([]*IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = NewIntVec(numCols, baseRing)
	}
	return &IntMatrix{numRows, numCols, rows, nil, nil, baseRing}
}

// NewIdIntMatrix returns an n by n identity matrix.
func NewIdIntMatrix(n int, baseRing *ring.Ring) *IntMatrix {
	m := NewIntMatrix(n, n, baseRing)
	for i, r := range m.rows {
		r.SetForce(i, 1)
	}
	return m
}

func NewIntMatrixFromRows(rows []*IntVec, baseRing *ring.Ring) *IntMatrix {
	return &IntMatrix{len(rows), rows[0].Size(), rows, nil, nil, baseRing}
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

func (m *IntMatrix) AsIntMatrix() *IntMatrix {
	return m
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

// ColCopy returns a copy of the i-th col.
func (m *IntMatrix) ColCopy(i int) *IntVec {
	colVec := NewIntVec(m.Rows(), m.baseRing)
	for j, row := range m.rows {
		colVec.SetCoeff(j, row.GetCoeff(i))
	}
	return colVec
}

// SubsectionCopy returns a subsection (copied) of this matrix.
func (m *IntMatrix) SubsectionCopy(rowStart, rowEnd int, colStart, colEnd int) *IntMatrix {
	subMatrix := NewIntMatrix(rowEnd-rowStart, colEnd-colStart, m.baseRing)
	for i := 0; i < rowEnd-rowStart; i++ {
		subMatrix.SetRow(i, m.RowView(rowStart+i).Slice(NewSlice(colStart, colEnd)))
	}
	return subMatrix
}

// SplitCols splits a matrix with k * d (where d is the poly. degree) cols into k many matrices with d cols,
// such that each row is represented by a single polynomial.
func (m *IntMatrix) SplitCols() []*IntMatrix {
	if m.Cols()%m.baseRing.N != 0 {
		panic("cannot split over columns")
	}
	numMatrices := m.Cols() / m.baseRing.N
	matrixRows := make([][]*IntVec, numMatrices)
	for i := 0; i < len(matrixRows); i++ {
		matrixRows[i] = make([]*IntVec, len(m.rows))
		for j, r := range m.rows {
			matrixRows[i][j] = r.UnderlyingPolysAsIntVecs()[i]
		}
	}
	matrices := make([]*IntMatrix, numMatrices)
	for i := 0; i < len(matrices); i++ {
		matrices[i] = NewIntMatrixFromRows(matrixRows[i], m.baseRing)
	}
	return matrices
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

// ExtendRows concatenates the matrices vertically.
func (m *IntMatrix) ExtendRows(b ImmutIntMatrix) {
	if m.Cols() != b.Cols() {
		panic("IntMatrix.ExtendRows cannot extend, invalid col size")
	}
	m.rows = append(m.rows, b.RowsView()...)
	m.numRows += b.Rows()
}

// ExtendCols concatenates the matrices horizontally.
func (m *IntMatrix) ExtendCols(b ImmutIntMatrix) {
	if m.Rows() != b.Rows() {
		panic("IntMatrix.ExtendCols cannot extend, invalid row size")
	}
	for i := 0; i < m.Rows(); i++ {
		m.rows[i].Append(b.RowView(i))
	}
	m.numCols += b.Cols()
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
func (m *IntMatrix) Copy() ImmutIntMatrix {
	rows := make([]*IntVec, m.numRows)
	for i, row := range m.rows {
		rows[i] = row.Copy()
	}
	return &IntMatrix{m.numRows, m.numCols, rows, m.unrebasedRef, m.rebaseRing, m.baseRing}
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
func (m *IntMatrix) RebaseRowsLossless(newRing RingParams) ImmutIntMatrix {
	rebasedRows := make([]*IntVec, len(m.rows))
	logging.LogShortExecution(fmt.Sprintf("%s.RebaseRowsLossless", m.SizeString()), "rebasing", func() interface{} {
		for i, row := range m.rows {
			rebasedRows[i] = row.RebaseLossless(newRing.BaseRing)
		}
		return nil
	})
	mp := NewIntMatrixFromRows(rebasedRows, newRing.BaseRing)
	mp.rebaseRing = &newRing
	mp.unrebasedRef = m
	return mp
}

// Scale scales the matrix by the given amount.
func (m *IntMatrix) Scale(factor uint64) MutIntMatrix {
	for _, row := range m.rows {
		row.Scale(factor)
	}
	return m
}

func (m *IntMatrix) ScaleCoeff(factors Coeff) MutIntMatrix {
	for _, row := range m.rows {
		row.ScaleCoeff(factors)
	}
	return m
}

// Neg negates this matrix.
func (m *IntMatrix) Neg() MutIntMatrix {
	for _, row := range m.rows {
		row.Neg()
	}
	return m
}

func (m *IntMatrix) IsUnset() bool {
	for _, r := range m.rows {
		if !r.IsUnset() {
			return false
		}
	}
	return true
}

// Transposed returns the transposed version of this matrix.
func (m *IntMatrix) Transposed() ImmutIntMatrix {
	if m.unrebasedRef != nil && m.rebaseRing != nil {
		return m.unrebasedRef.Transposed().RebaseRowsLossless(*m.rebaseRing)
	}
	mt := NewIntMatrix(m.Cols(), m.Rows(), m.baseRing)
	if m.IsUnset() {
		return mt
	}
	procName := fmt.Sprintf("%s.Transposed", m.SizeString())
	e := logging.LogExecShortStart(procName, "transposing")
	defer e.LogExecEnd()
	for lvl := 0; lvl < len(m.baseRing.Modulus); lvl++ {
		// pull in the level
		A := make([][]uint64, m.Rows())
		for i, row := range m.rows {
			A[i] = row.GetWholeLevel(lvl)
		}
		// transpose the level
		At := make([][]uint64, m.Cols())
		for i := 0; i < m.Cols(); i++ {
			At[i] = make([]uint64, m.Rows())
			for j := 0; j < m.Rows(); j++ {
				At[i][j] = A[j][i]
			}
		}
		// move out
		for i, r := range At {
			mt.RowView(i).SetWholeLevel(lvl, r)
		}
	}
	return mt
}

// Hadamard performs coefficient-wise multiplication.
func (m *IntMatrix) Hadamard(b ImmutIntMatrix) ImmutIntMatrix {
	for i, r := range m.rows {
		r.Hadamard(b.RowView(i))
	}
	return m
}

// Max returns the largest element of the matrix.
func (m *IntMatrix) Max(q *big.Int) *big.Int {
	max := m.rows[0].Max(q)
	for _, r := range m.rows {
		c := r.Max(q)
		if c.Cmp(max) > 0 {
			max = c
		}
	}
	return max
}

// Eq returns true iff two matrices are equal in their elements.
func (m *IntMatrix) Eq(b ImmutIntMatrix) bool {
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

// MulVecSlow performs a matrix-vector multiplication using an unoptimized multiplication procedure. Use `MulVec` instead!
func (m *IntMatrix) MulVecSlow(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	if m.IsUnset() {
		return NewIntVec(m.Rows(), m.baseRing)
	}
	return logging.LogShortExecution(fmt.Sprintf("%s.MulVec", m.SizeString()), "multiplying", func() interface{} {
		if m.unrebasedRef != nil && v.unrebasedRef != nil {
			return m.unrebasedRef.MulVec(v.unrebasedRef)
		}
		out := NewIntVec(m.Rows(), m.baseRing)
		for i, row := range m.rows {
			dotResult := row.Dot(v)
			out.SetCoeff(i, dotResult)
		}
		return out
	}).(*IntVec)
}

// MulVec performs a matrix-vector multiplication.
func (m *IntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	if m.IsUnset() {
		return NewIntVec(m.Cols(), m.baseRing)
	}
	procName := fmt.Sprintf("%s.MulVec", m.SizeString())
	return logging.LogExecution(procName, "multiplying",
		func() interface{} {
			if m.unrebasedRef != nil && v.unrebasedRef != nil {
				logging.Log(procName, "using unrebased references directly")
				res := MulVecBlazingFast(m.unrebasedRef, v.unrebasedRef, m.unrebasedRef.baseRing)
				// rebase back
				return res.RebaseLossless(v.baseRing)
			} else if m.unrebasedRef != nil && v.unrebasedRef == nil && v.Size() == m.unrebasedRef.baseRing.N {
				logging.Log(procName, "unrebasing the vector")
				v.unrebasedRef = v.MergedPolys(m.unrebasedRef.baseRing)
				res := MulVecBlazingFast(m.unrebasedRef, v.unrebasedRef, m.unrebasedRef.baseRing)
				// rebase back
				return res.RebaseLossless(v.baseRing)
			}
			logging.Log(procName, "unrebase not possible")
			return MulVecBlazingFast(m, v, m.baseRing)
		}).(*IntVec)
}

// MulVecTranspose performs a matrix-vector multiplication with the transpose of this matrix.
func (m *IntMatrix) MulVecTranspose(v *IntVec) *IntVec {
	if m.Rows() != v.Size() {
		panic("IntMatrix.MulVecTranspose sizes incorrect")
	}
	if m.IsUnset() {
		return NewIntVec(m.Cols(), m.baseRing)
	}
	procName := fmt.Sprintf("%s.MulVecTranspose", m.SizeString())
	return logging.LogExecution(procName, "multiplying",
		func() interface{} {
			if m.unrebasedRef != nil && v.unrebasedRef != nil {
				logging.Log(procName, "using unrebased references directly")
				res := MulVecTransposeBlazingFast(m.unrebasedRef, v.unrebasedRef, m.unrebasedRef.baseRing)
				// rebase back
				return res.RebaseLossless(v.baseRing)
			} else if m.unrebasedRef != nil && v.unrebasedRef == nil && v.Size() == m.unrebasedRef.baseRing.N {
				logging.Log(procName, "unrebasing the vector")
				v.unrebasedRef = v.MergedPolys(m.unrebasedRef.baseRing)
				res := MulVecTransposeBlazingFast(m.unrebasedRef, v.unrebasedRef, m.unrebasedRef.baseRing)
				// rebase back
				return res.RebaseLossless(v.baseRing)
			}
			logging.Log(procName, "unrebase not possible")
			return MulVecTransposeBlazingFast(m, v, m.baseRing)
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

func (m *IntMatrix) Cleanup() {
	for _, r := range m.rows {
		r.Cleanup()
	}
}
