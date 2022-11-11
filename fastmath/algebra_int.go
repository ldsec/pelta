package fastmath

import (
	"fmt"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

type IntVec struct {
	size     int
	polys    []Poly
	mod      *big.Int
	baseRing *ring.Ring
}

func NewIntVec(size int, baseRing *ring.Ring) *IntVec {
	numPolys := int(size/baseRing.N) + 1
	if size%baseRing.N == 0 {
		numPolys -= 1
	}
	polys := make([]Poly, numPolys)
	for i := 0; i < len(polys); i++ {
		polys[i] = *NewZeroPoly(baseRing)
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

func NewIntVecFromPolys(polys []Poly, size int, baseRing *ring.Ring) *IntVec {
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
func (v *IntVec) UnderlyingPolys() []Poly {
	return v.polys
}

// SetUnderlyingPolys can be used to update the polynomials that are used to represent this vector.
func (v *IntVec) SetUnderlyingPolys(polys []Poly) {
	v.polys = polys
}

// Set updates the given element of this vector.
func (v *IntVec) Set(index int, newValue uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].Set(coeffIndex, newValue)
}

// Scale scales the vector by the given amount.
func (v *IntVec) Scale(factor uint64) *IntVec {
	for _, p := range v.polys {
		p.Scale(factor)
	}
	return v
}

// Neg negates the polynomial.
func (v *IntVec) Neg() *IntVec {
	for _, p := range v.polys {
		p.Neg()
	}
	return v
}

// Add adds two integer vectors.
func (v *IntVec) Add(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Add sizes do not match")
	}
	for i, p := range v.polys {
		p.Add(&r.polys[i])
	}
	return v
}

// Dot returns the dot product of the given two vectors.
func (v *IntVec) Dot(r *IntVec) uint64 {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	preSum := NewZeroPoly(v.baseRing)
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		v.baseRing.MulCoeffsAndAdd(a.ref, b.ref, preSum.ref)
	}
	return preSum.SumCoeffsLimited(0, v.size, v.mod)
}

// MulAddElems multiplies the elements and adds it to the coefficients of the
// given `out` polynomial.
func (v *IntVec) MulAddElems(r *IntVec, out *Poly) {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		v.baseRing.MulCoeffsAndAdd(a.ref, b.ref, out.ref)
	}
}

// Eq checks the equality between two integer vectors.
func (v *IntVec) Eq(r *IntVec) bool {
	if v.Size() != r.Size() {
		return false
	}
	// The underlying polynomial # might change even when the sizes are equal.
	// Take the minimum.
	minPolySize := len(v.polys)
	if len(r.polys) < minPolySize {
		minPolySize = len(r.polys)
	}
	for i := 0; i < minPolySize; i++ {
		if !v.polys[i].EqLevel(0, &r.polys[i]) {
			return false
		}
	}
	return true
}

// Copy copies the vector.
func (v *IntVec) Copy() *IntVec {
	polys := make([]Poly, len(v.polys))
	for i, p := range v.polys {
		polys[i] = *p.Copy()
	}
	return &IntVec{v.size, polys, v.mod, v.baseRing}
}

// String returns the string representation of this integer vector.
func (v *IntVec) String() string {
	s := fmt.Sprintf("IntVec[%d]{", v.Size())
	elemStrs := make([]string, 0, v.Size())
	for i := 0; i < len(v.polys); i++ {
		for j := 0; j < v.polys[i].N(); j++ {
			elemStrs = append(elemStrs, fmt.Sprintf("%d", v.polys[i].Get(j, 0)))
		}
	}
	return s + strings.Join(elemStrs[:v.Size()], ",") + ",...}"
}

// Append appends the contents of the given vector into this one.
func (v *IntVec) Append(r *IntVec) *IntVec {
	v.size = v.size + r.size
	v.polys = append(v.polys, r.polys...)
	return v
}

// Reduce reduces the elements of this vector by the given mod.
func (v *IntVec) Reduce(mod *big.Int) *IntVec {
	for _, p := range v.polys {
		for i := 0; i < p.N(); i++ {
			reducedVal := big.NewInt(0).Mod(big.NewInt(int64(p.Get(i, 0))), mod)
			p.Set(i, reducedVal.Uint64())
		}
	}
	return v
}

type IntMatrix struct {
	numRows  int
	numCols  int
	rows     []IntVec
	mod      *big.Int
	baseRing *ring.Ring
}

func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) *IntMatrix {
	rows := make([]IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = *NewIntVec(numCols, baseRing)
	}
	return &IntMatrix{numRows, numCols, rows, baseRing.ModulusAtLevel[0], baseRing}
}

func NewIntMatrixFromSlice(elems [][]uint64, baseRing *ring.Ring) *IntMatrix {
	numRows := len(elems)
	numCols := len(elems[0])
	m := NewIntMatrix(numRows, numCols, baseRing)
	m.PopulateRows(func(i int) IntVec {
		return *NewIntVecFromSlice(elems[i], baseRing)
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

// RowView returns a refernce to the i-th row.
func (m *IntMatrix) RowView(i int) *IntVec {
	return &m.rows[i]
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
func (m *IntMatrix) SetRow(row int, newRow IntVec) {
	if row >= m.Rows() {
		panic("IntMatrix.SetRows index incorrect")
	}
	m.rows[row] = newRow
}

// AppendRow appends a row into the vector.
func (m *IntMatrix) AppendRow(v IntVec) {
	if m.Cols() != v.Size() {
		panic("IntMatrix.AppendRow cannot append row, invalid size")
	}
	m.rows = append(m.rows, v)
	m.numRows += 1
}

// ExtendRows concatenates the matrices vertically.
func (m *IntMatrix) ExtendRows(b *IntMatrix) *IntMatrix {
	if m.Cols() != b.Cols() {
		panic("IntMatrix.ExtendRows cannot extend, invalid size")
	}
	m.rows = append(m.rows, b.rows...)
	m.numRows += b.Rows()
	return m
}

// ExtendCols concatentes the matrices horizontally.
func (m *IntMatrix) ExtendCols(b *IntMatrix) *IntMatrix {
	if m.Rows() != b.Rows() {
		panic("IntMatrix.ExtendCols cannot extend, invalid size")
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
func (m *IntMatrix) PopulateRows(f func(int) IntVec) {
	for i := 0; i < m.Rows(); i++ {
		m.rows[i] = f(i)
	}
}

// Transposed returns the transposed version of this matrix.
func (m *IntMatrix) Transposed() *IntMatrix {
	newRows := make([]IntVec, m.Cols())
	for i := 0; i < len(newRows); i++ {
		newRows[i] = *m.ColCopy(i)
	}
	return &IntMatrix{m.numCols, m.numRows, newRows, m.mod, m.baseRing}
}

// MulVec performs a matrix-vector multiplication.
func (m *IntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	out := NewIntVec(m.Rows(), m.baseRing)
	dotResult := NewZeroPoly(m.baseRing)
	for i, row := range m.rows {
		row.MulAddElems(v, dotResult)
		out.Set(i, dotResult.SumCoeffsLimited(0, row.Size(), m.mod))
		dotResult.ref.Zero()
	}
	return out
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
			out.Set(i, j, p.Dot(q))
		}
	}
	return out
}

// Scale scales the matrix by the given amount.
func (m *IntMatrix) Scale(factor uint64) {
	for _, row := range m.rows {
		row.Scale(factor)
	}
}

// Copy returns a copy of this matrix.
func (m *IntMatrix) Copy() *IntMatrix {
	rows := make([]IntVec, m.numRows)
	for i, row := range m.rows {
		rows[i] = *row.Copy()
	}
	return &IntMatrix{m.numRows, m.numCols, rows, m.mod, m.baseRing}
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

// String returns a string representation of the matrix.
func (m *IntMatrix) String() string {
	s := fmt.Sprintf("IntMatrix[%d,%d]{\n", m.Rows(), m.Cols())
	for _, row := range m.rows {
		s += "\t" + row.String() + "\n"
	}
	return s + ", ...}"
}
