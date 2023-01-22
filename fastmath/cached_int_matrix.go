package fastmath

import (
	"fmt"
	"strings"

	"github.com/ldsec/codeBase/commitment/logging"
)

// CachedIntMatrix represents a mutable int matrix whose transpose is generated on the fly as
// the modifications are performed on the original.
type CachedIntMatrix struct {
	A  MutIntMatrix
	At MutIntMatrix
}

func NewCachedIntMatrix(m MutIntMatrix) *CachedIntMatrix {
	return &CachedIntMatrix{m, m.Transposed().(MutIntMatrix)}
}

// Rows returns the number of rows.
func (m *CachedIntMatrix) Rows() int {
	return m.A.Rows()
}

// Cols returns the number of cols.
func (m *CachedIntMatrix) Cols() int {
	return m.A.Cols()
}

// RowsView returns a list of references to the rows of this matrix.
func (m *CachedIntMatrix) RowsView() []*IntVec {
	return m.A.RowsView()
}

// RowView returns a reference to the i-th row.
func (m *CachedIntMatrix) RowView(i int) *IntVec {
	return m.A.RowView(i)
}

// ColCopy returns a copy of the i-th col.
func (m *CachedIntMatrix) ColView(i int) *IntVec {
	return m.At.RowView(i)
}

// SubsectionCopy returns a subsection (copied) of this matrix.
func (m *CachedIntMatrix) SubsectionCopy(rowStart, rowEnd int, colStart, colEnd int) *IntMatrix {
	return m.A.SubsectionCopy(rowStart, rowEnd, colStart, colEnd)
}

// GetCoeff returns the element at the given coordinates.
func (m *CachedIntMatrix) GetCoeff(row, col int) Coeff {
	return m.A.GetCoeff(row, col)
}

// Get returns the element at the given coordinates.
func (m *CachedIntMatrix) GetLevel(row, col, level int) uint64 {
	return m.A.GetLevel(row, col, level)
}

// SetForce updates the given element of this matrix.
func (m *CachedIntMatrix) SetForce(row, col int, newValue uint64) {
	m.A.SetForce(row, col, newValue)
	m.At.SetForce(col, row, newValue)
}

// SetCoeff updates the given element of this matrix.
func (m *CachedIntMatrix) SetCoeff(row, col int, newCoeff Coeff) {
	m.A.SetCoeff(row, col, newCoeff)
	m.At.SetCoeff(col, row, newCoeff)
}

// SetRow updates the given row of this matrix.
func (m *CachedIntMatrix) SetRow(row int, newRow *IntVec) {
	m.A.SetRow(row, newRow)
	m.At.SetCol(row, newRow)
}

// SetCol updates the given col of this matrix.
func (m *CachedIntMatrix) SetCol(col int, newCol *IntVec) {
	m.A.SetCol(col, newCol)
	m.At.SetRow(col, newCol)
}

// ExtendRows concatenates the matrices vertically.
func (m *CachedIntMatrix) ExtendRows(b ImmutIntMatrix) {
	m.A.ExtendRows(b)
	m.At.ExtendCols(b.Transposed())
}

// ExtendCols concatenates the matrices horizontally.
func (m *CachedIntMatrix) ExtendCols(b ImmutIntMatrix) {
	m.A.ExtendCols(b)
	m.At.ExtendRows(b.Transposed())
}

// Populate is used to initialize the elements of this matrix.
func (m *CachedIntMatrix) Populate(f func(int, int) uint64) {
	m.A.Populate(f)
	m.At.Populate(func(i, j int) uint64 { return f(j, i) })
}

// Copy returns a copy of this matrix.
func (m *CachedIntMatrix) Copy() ImmutIntMatrix {
	return &CachedIntMatrix{m.A.Copy().(*IntMatrix), m.At.Copy().(*IntMatrix)}
}

// String returns a string representation of the matrix.
func (m *CachedIntMatrix) String() string {
	return strings.Replace(m.A.String(), "IntMatrix", "CachedIntMatrix", 1)
}

// SizeString returns a string representation of the matrix's dimensions.
func (m *CachedIntMatrix) SizeString() string {
	return strings.Replace(m.A.SizeString(), "IntMatrix", "CachedIntMatrix", 1)
}

// RebaseRowsLossless rebases each row.
func (m *CachedIntMatrix) RebaseRowsLossless(newRing RingParams) ImmutIntMatrix {
	return &CachedIntMatrix{
		m.A.RebaseRowsLossless(newRing).(MutIntMatrix),
		m.At.RebaseRowsLossless(newRing).(MutIntMatrix),
	}
}

// Scale scales the matrix by the given amount.
func (m *CachedIntMatrix) Scale(factor uint64) MutIntMatrix {
	m.A.Scale(factor)
	m.At.Scale(factor)
	return m
}

func (m *CachedIntMatrix) ScaleCoeff(factors Coeff) MutIntMatrix {
	m.A.ScaleCoeff(factors)
	m.At.ScaleCoeff(factors)
	return m
}

// Neg negates this matrix.
func (m *CachedIntMatrix) Neg() MutIntMatrix {
	m.A.Neg()
	m.At.Neg()
	return m
}

// Transposed returns the transposed version of this matrix.
func (m *CachedIntMatrix) Transposed() ImmutIntMatrix {
	procName := fmt.Sprintf("%s.Transposed (cached)", m.SizeString())
	e := logging.LogExecShortStart(procName, "transposing")
	defer e.LogExecEnd()
	return m.At
}

// Hadamard performs coefficient-wise multiplication.
func (m *CachedIntMatrix) Hadamard(b ImmutIntMatrix) ImmutIntMatrix {
	m.A.Hadamard(b)
	m.At.Hadamard(b.Transposed())
	return m
}

// Max returns the largest element of the matrix.
func (m *CachedIntMatrix) Max() uint64 {
	return m.A.Max()
}

// Eq returns true iff two matrices are equal in their elements.
func (m *CachedIntMatrix) Eq(b ImmutIntMatrix) bool {
	return m.A.Eq(b)
}

// MulVec performs a matrix-vector multiplication.
func (m *CachedIntMatrix) MulVec(v *IntVec) *IntVec {
	return m.A.MulVec(v)
}

func (m *CachedIntMatrix) AsIntMatrix() *IntMatrix {
	return m.A.AsIntMatrix()
}

func (m *CachedIntMatrix) Cleanup() {
	m.A.Cleanup()
	m.At.Cleanup()
}
