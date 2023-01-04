package fastmath

import (
	"fmt"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

type IntVec struct {
	size     int
	polys    []*Poly
	mod      *big.Int
	baseRing *ring.Ring
}

func NewIntVec(size int, baseRing *ring.Ring) *IntVec {
	numPolys := int(size/baseRing.N) + 1
	if size%baseRing.N == 0 {
		numPolys -= 1
	}
	polys := make([]*Poly, numPolys)
	for i := 0; i < len(polys); i++ {
		polys[i] = NewPoly(baseRing)
	}
	return &IntVec{size, polys, baseRing.ModulusAtLevel[0], baseRing}
}

func NewIntVecFromSlice(slice []uint64, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(len(slice), baseRing)
	for i := 0; i < len(slice); i++ {
		v.Set(i, slice[i])
	}
	return v
}

func NewIntVecFromPolys(polys []*Poly, size int, baseRing *ring.Ring) *IntVec {
	return &IntVec{size, polys, baseRing.ModulusAtLevel[0], baseRing}
}

// Size returns the size of this vector.
func (v *IntVec) Size() int {
	return v.size
}

// Populate is used to initialize the elements of this vector.
func (v *IntVec) Populate(f func(int) uint64) {
	for i := 0; i < v.size; i++ {
		val := f(i)
		v.Set(i, val)
	}
}

// Get returns the element at the given index.
func (v *IntVec) Get(index int) uint64 {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	return v.polys[polyIndex].Get(coeffIndex, 0)
}

// UnderlyingPolys returns the polynomials that are being used to represent this vector.
func (v *IntVec) UnderlyingPolys() []*Poly {
	return v.polys
}

func (v *IntVec) UnderlyingPolysAsPolyVec() *PolyVec {
	pV := NewPolyVec(len(v.polys), v.BaseRing())
	for i, p := range v.polys {
		pV.Set(i, p)
	}
	return pV
}

func (v *IntVec) UnderlyingPolysAsPolyNTTVec() *PolyNTTVec {
	pV := NewPolyVec(len(v.polys), v.BaseRing()).NTT()
	for i, p := range v.polys {
		pV.Set(i, ForceNTT(p))
	}
	return pV
}

// SetUnderlyingPolys can be used to update the polynomials that are used to represent this vector.
func (v *IntVec) SetUnderlyingPolys(polys []*Poly) {
	v.polys = polys
}

// Set updates the given element of this vector.
func (v *IntVec) Set(index int, newValue uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].Set(coeffIndex, newValue)
}

// Copy copies the vector.
func (v *IntVec) Copy() *IntVec {
	polys := make([]*Poly, len(v.polys))
	for i, p := range v.polys {
		polys[i] = p.Copy()
	}
	return &IntVec{v.size, polys, v.mod, v.baseRing}
}

// String returns the string representation of this integer vector.
func (v *IntVec) String() string {
	s := fmt.Sprintf("IntVec[%d]{", v.Size())
	elemStrs := make([]string, 0, v.Size())
	for i := 0; i < v.Size(); i++ {
		elemStrs = append(elemStrs, fmt.Sprintf("%d", v.Get(i)))
	}
	return s + strings.Join(elemStrs, ",") + "}"
}

// Append appends the contents of the given vector into this one.
func (v *IntVec) Append(r *IntVec) *IntVec {
	// Update the size.
	newSize := v.size + r.size
	// Optimization.
	if v.Size()%v.baseRing.N == 0 {
		v.size = newSize
		v.polys = append(v.polys, r.polys...)
		return v
	}
	// Extend the number of underlying polynomials.
	numPolys := int(newSize/v.baseRing.N) + 1
	if newSize%v.baseRing.N == 0 {
		numPolys -= 1
	}
	oldSize := v.size
	v.size = newSize
	for i := len(v.polys); i < numPolys; i++ {
		v.polys = append(v.polys, NewPoly(v.baseRing))
	}
	// Move the elements in.
	for i := 0; i < r.Size(); i++ {
		v.Set(oldSize+i, r.Get(i))
	}
	return v
}

// SliceCopy returns a (copied) slice of this vector.
func (v *IntVec) SliceCopy(start, end int) *IntVec {
	subVec := NewIntVec(end-start, v.baseRing)
	for i := 0; i < end-start; i++ {
		subVec.Set(i, v.Get(start+i))
	}
	return subVec
}

// RebaseLossless rebases every underlying polynomial, splitting them when necessary so that
// no coefficient will be discarded.
func (v *IntVec) RebaseLossless(newRing RingParams, level int) *IntVec {
	if v.baseRing.N <= newRing.D || v.baseRing.N%newRing.D != 0 {
		panic(fmt.Sprintf("cannot rebase lossless %d to %d", v.baseRing.N, newRing.D))
	}
	// Clean up the redundant polynomials.
	// for len(v.polys) > 0 && v.polys[len(v.polys)-1].IsZero() {
	// 	v.polys = v.polys[0 : len(v.polys)-1]
	// }
	// Each underlying polynomial will be represented by this many rebased polynomials.
	splitsPerPoly := v.baseRing.N / newRing.D
	newPolys := []*Poly{}
	for _, p := range v.polys {
		for i := 0; i < splitsPerPoly; i++ {
			// Extract d (i.e., degree of new base ring) many coefficients.
			newCoeffs := p.ref.Coeffs[level][i*newRing.D : (i+1)*newRing.D]
			newPoly := NewPoly(newRing.BaseRing)
			for j := 0; j < newPoly.N(); j++ {
				newPoly.Set(j, newCoeffs[j])
			}
			newPolys = append(newPolys, newPoly)
		}
	}
	// Set the underlying polynomials.
	v.SetUnderlyingPolys(newPolys)
	// Update the base ring pointer.
	v.baseRing = newRing.BaseRing
	v.mod = newRing.BaseRing.ModulusAtLevel[0]
	return v
}

func (v *IntVec) All(pred func(el uint64) bool) bool {
	for i := 0; i < v.Size(); i++ {
		if !pred(v.Get(i)) {
			return false
		}
	}
	return true
}

// BaseRing returns the polynomial ring over which this integer vector is defined.
func (v *IntVec) BaseRing() *ring.Ring {
	return v.baseRing
}

type IntMatrix struct {
	numRows  int
	numCols  int
	rows     []*IntVec
	mod      *big.Int
	baseRing *ring.Ring
}

func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) *IntMatrix {
	rows := make([]*IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = NewIntVec(numCols, baseRing)
	}
	return &IntMatrix{numRows, numCols, rows, baseRing.ModulusAtLevel[0], baseRing}
}

// NewIdIntMatrix returns an n by n identity matrix.
func NewIdIntMatrix(numRows int, baseRing *ring.Ring) *IntMatrix {
	m := NewIntMatrix(numRows, numRows, baseRing)
	for i, r := range m.RowsView() {
		r.Set(i, 1)
	}
	return m
}

func NewIntMatrixFromRows(rows []*IntVec, baseRing *ring.Ring) *IntMatrix {
	return &IntMatrix{len(rows), rows[0].Size(), rows, baseRing.ModulusAtLevel[0], baseRing}
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
		colVec.Set(j, row.Get(i))
	}
	return colVec
}

// Get returns the element at the given coordinates.
func (m *IntMatrix) Get(row, col int) uint64 {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Get indices incorrect")
	}
	return m.rows[row].Get(col)
}

// Set updates the given element of this matrix.
func (m *IntMatrix) Set(row, col int, newValue uint64) {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	m.rows[row].Set(col, newValue)
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
		row.Set(col, newCol.Get(i))
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
	return &IntMatrix{m.numRows, m.numCols, rows, m.mod, m.baseRing}
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
	s := fmt.Sprintf("IntMatrix[%d,%d]\n", m.Rows(), m.Cols())
	return s
}

// RebaseRowsLossless rebases each row.
func (m *IntMatrix) RebaseRowsLossless(newRing RingParams, level int) *IntMatrix {
	m.baseRing = newRing.BaseRing
	for _, row := range m.rows {
		row.RebaseLossless(newRing, level)
	}
	return m
}
