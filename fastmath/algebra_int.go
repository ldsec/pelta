package fastmath

import (
	"fmt"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

type IntVec struct {
	size     int
	polys    []Poly
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
	return &IntVec{size, polys, baseRing}
}

func NewIntVecFromSlice(slice []uint64, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(len(slice), baseRing)
	for i := 0; i < len(slice); i++ {
		v.Set(i, slice[i])
	}
	return v
}

func NewIntVecFromPolys(polys []Poly, size int, baseRing *ring.Ring) *IntVec {
	return &IntVec{size, polys, baseRing}
}

func (v *IntVec) Size() int {
	return v.size
}

func (v *IntVec) Populate(f func(int) uint64) {
	for i := 0; i < v.size; i++ {
		val := f(i)
		v.Set(i, val)
	}
}

func (v *IntVec) Get(index int) uint64 {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	return v.polys[polyIndex].Get(coeffIndex, 0)
}

// UnderlyingPolys returns the polynomials that are being used to represent this vector.
func (v *IntVec) UnderlyingPolys() []Poly {
	return v.polys
}

func (v *IntVec) SetUnderlyingPolys(polys []Poly) {
	v.polys = polys
}

func (v *IntVec) Set(index int, newValue uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].Set(coeffIndex, newValue)
}

// Scale scales the vector by the given amount.
func (v *IntVec) Scale(factor uint64) {
	for _, p := range v.polys {
		p.Scale(factor)
	}
}

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
	sum := uint64(0)
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		c := NewZeroPoly(v.baseRing)
		v.baseRing.MulCoeffsAndAdd(a.ref, b.ref, c.ref)
		sum += c.SumCoeffs(0)
	}
	return sum
}

// Eq checks the equality between two integer vectors.
func (v *IntVec) Eq(r *IntVec) bool {
	if v.size != r.size {
		return false
	}
	for i := 0; i < len(v.polys); i++ {
		if !v.polys[i].EqLevel(0, &r.polys[i]) {
			return false
		}
	}
	return true
}

func (v *IntVec) Copy() *IntVec {
	polys := make([]Poly, len(v.polys))
	for i, p := range v.polys {
		polys[i] = *p.Copy()
	}
	return &IntVec{size: v.size, polys: polys, baseRing: v.baseRing}
}

// String returns the string representation of this integer vector.
func (v *IntVec) String() string {
	s := fmt.Sprintf("IntVec[%d]{", v.Size())
	elemStrs := make([]string, 0, v.size)
	for i := 0; i < len(v.polys); i++ {
		for j := 0; j < v.polys[i].N(); j++ {
			elemStrs = append(elemStrs, fmt.Sprintf("%d", v.polys[i].Get(j, 0)))
		}
	}
	return s + strings.Join(elemStrs[:v.size], ",") + ",...}"
}

type IntMatrix struct {
	numRows  int
	numCols  int
	rows     []IntVec
	baseRing *ring.Ring
}

func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) *IntMatrix {
	rows := make([]IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = *NewIntVec(numCols, baseRing)
	}
	return &IntMatrix{numRows, numCols, rows, baseRing}
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

func (m *IntMatrix) Rows() int {
	return m.numRows
}

func (m *IntMatrix) Cols() int {
	return m.numCols
}

func (m *IntMatrix) RowView(i int) *IntVec {
	return &m.rows[i]
}

func (m *IntMatrix) ColCopy(i int) *IntVec {
	colVec := NewIntVec(m.Rows(), m.baseRing)
	for j, row := range m.rows {
		colVec.Set(j, row.Get(i))
	}
	return colVec
}

func (m *IntMatrix) Get(row, col int) uint64 {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Get indices incorrect")
	}
	return m.rows[row].Get(col)
}

func (m *IntMatrix) Set(row, col int, newValue uint64) {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	m.rows[row].Set(col, newValue)
}

func (m *IntMatrix) SetRow(row int, newRow IntVec) {
	if row >= m.Rows() {
		panic("IntMatrix.SetRows index incorrect")
	}
	m.rows[row] = newRow
}

func (m *IntMatrix) Populate(f func(int, int) uint64) {
	for i, row := range m.rows {
		row.Populate(
			func(j int) uint64 {
				return f(i, j)
			})
	}
}

func (m *IntMatrix) PopulateRows(f func(int) IntVec) {
	for i := 0; i < m.Rows(); i++ {
		m.rows[i] = f(i)
	}
}

func (m *IntMatrix) Transposed() *IntMatrix {
	newRows := make([]IntVec, m.Cols())
	for i := 0; i < len(newRows); i++ {
		newRows[i] = *m.ColCopy(i)
	}
	return &IntMatrix{m.numCols, m.numRows, newRows, m.baseRing}
}

func (m *IntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	out := NewIntVec(m.Rows(), m.baseRing)
	for i, row := range m.rows {
		out.Set(i, row.Dot(v))
	}
	return out
}

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

func (m *IntMatrix) Copy() *IntMatrix {
	rows := make([]IntVec, m.numRows)
	for i, row := range m.rows {
		rows[i] = *row.Copy()
	}
	return &IntMatrix{m.numRows, m.numCols, rows, m.baseRing}
}

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

func (m *IntMatrix) String() string {
	s := fmt.Sprintf("IntMatrix[%d,%d]{\n", m.Rows(), m.Cols())
	for _, row := range m.rows {
		s += "\t" + row.String() + "\n"
	}
	return s + ", ...}"
}