package fastmath

import (
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
)

type Poly struct {
	ref      *ring.Poly
	baseRing *ring.Ring
}

// NewZeroPoly returns a zero polynomial.
func NewZeroPoly(baseRing *ring.Ring) Poly {
	return Poly{baseRing.NewPoly(), baseRing}
}

// NewOnePoly returns a one polynomial scaled with the given factor.
func NewOnePoly(scale uint64, baseRing *ring.Ring) Poly {
	p := NewZeroPoly(baseRing)
	p.Set(0, scale)
	return p
}

// Set sets the coefficient of this polynomial at every level to the given value.
func (p *Poly) Set(index int, value uint64) {
	for level := 0; level < len(p.ref.Coeffs); level++ {
		p.ref.Coeffs[level][index] = value
	}
}

// Get returns the coefficient of this polynomial.
func (p *Poly) Get(index int, level int) uint64 {
	return p.ref.Coeffs[level][index]
}

// SumCoeffs returns the sum of the coefficients of this polynomial.
func (p *Poly) SumCoeffs(level int) uint64 {
	logN := int(math.Log2(float64(p.baseRing.N)))
	tmp := p.Copy()
	for i := 0; i < logN; i++ {
		tmp2 := NewZeroPoly(p.baseRing)
		p.baseRing.Shift(tmp.ref, 1<<i, tmp2.ref)
		p.baseRing.Add(tmp.ref, tmp2.ref, tmp.ref)
	}
	return tmp.ref.Coeffs[level][0]
}

// Scale scales this polynomial with the given scalar factor.
func (p *Poly) Scale(factor uint64) {
	p.baseRing.MulScalar(p.ref, factor, p.ref)
}

// NTT converts this polynomial to its NTT domain.
func (p *Poly) NTT() {
	p.baseRing.NTT(p.ref, p.ref)
}

// InvNTT converts this polynomial back to its poly domain.
func (p *Poly) InvNTT() {
	p.baseRing.InvNTT(p.ref, p.ref)
}

// Add adds two polynomials.
func (p *Poly) Add(other *Poly) {
	p.baseRing.Add(p.ref, other.ref, p.ref)
}

// MulCoeffs multiplies the coefficients of two polynomials.
func (p *Poly) MulCoeffs(other *Poly) {
	p.baseRing.MulCoeffs(p.ref, other.ref, p.ref)
}

// Pow takes the `exp`-th power of the coefficients modulo `mod`.
func (p *Poly) PowModCoeffs(exp uint64, mod uint64) {
	for i := 0; i < p.baseRing.N; i++ {
		newCoeff := ring.ModExp(p.Get(i, 0), exp, mod)
		p.Set(i, newCoeff)
	}
}

// Copy returns a copy of this polynomial.
func (p *Poly) Copy() Poly {
	return Poly{p.ref.CopyNew(), p.baseRing}
}

type PolyVec struct {
	elems    []Poly
	baseRing *ring.Ring
}

func NewPolyVec(size int, baseRing *ring.Ring) PolyVec {
	elems := make([]Poly, size)
	for i := 0; i < size; i++ {
		elems[i] = NewZeroPoly(baseRing)
	}
	return PolyVec{elems, baseRing}
}

type PolyMatrix struct {
	elems    [][]Poly
	baseRing *ring.Ring
}

func NewPolyMatrix(rows, cols int, baseRing *ring.Ring) PolyMatrix {
	elems := make([][]Poly, rows)
	for i := 0; i < rows; i++ {
		elems[i] = make([]Poly, cols)
		for j := 0; j < cols; j++ {
			elems[i][j] = NewZeroPoly(baseRing)
		}
	}
	return PolyMatrix{elems, baseRing}
}

type IntVec struct {
	size     int
	polys    []Poly
	baseRing *ring.Ring
}

func NewIntVec(size int, baseRing *ring.Ring) IntVec {
	numPolys := int(size/baseRing.N) + 1
	polys := make([]Poly, numPolys)
	for i := 0; i < len(polys); i++ {
		polys[i] = NewZeroPoly(baseRing)
	}
	return IntVec{size, polys, baseRing}
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

func (v *IntVec) Set(index int, newValue uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].Set(coeffIndex, newValue)
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

type IntMatrix struct {
	numRows  int
	numCols  int
	rows     []IntVec
	baseRing *ring.Ring
}

func NewIntMatrix(numRows, numCols int, baseRing *ring.Ring) IntMatrix {
	rows := make([]IntVec, numRows)
	for i := 0; i < len(rows); i++ {
		rows[i] = NewIntVec(numCols, baseRing)
	}
	return IntMatrix{numRows, numCols, rows, baseRing}
}

func (m *IntMatrix) Rows() int {
	return m.numRows
}

func (m *IntMatrix) Cols() int {
	return m.numCols
}

func (m *IntMatrix) Row(i int) IntVec {
	return m.rows[i]
}

func (m *IntMatrix) Col(i int) IntVec {
	colVec := NewIntVec(m.Rows(), m.baseRing)
	for j, row := range m.rows {
		colVec.Set(j, row.Get(i))
	}
	return colVec
}

func (m *IntMatrix) Get(row, col int) uint64 {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	return m.rows[row].Get(col)
}

func (m *IntMatrix) Set(row, col int, newValue uint64) {
	if row >= m.Rows() || col >= m.Cols() {
		panic("IntMatrix.Set indices incorrect")
	}
	m.rows[row].Set(col, newValue)
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

func (m *IntMatrix) Transposed() IntMatrix {
	newRows := make([]IntVec, m.Cols())
	for i := 0; i < len(newRows); i++ {
		newRows[i] = m.Col(i)
	}
	return IntMatrix{m.numCols, m.numRows, newRows, m.baseRing}
}

func (m *IntMatrix) MulVec(v *IntVec) IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	out := NewIntVec(m.Rows(), m.baseRing)
	for i, row := range m.rows {
		out.Set(i, row.Dot(v))
	}
	return out
}

func (m *IntMatrix) MulMat(b *IntMatrix) IntMatrix {
	if m.Cols() != b.Rows() {
		panic("IntMatrix.MulMat sizes incorrect")
	}
	out := NewIntMatrix(m.Rows(), b.Cols(), m.baseRing)
	for i := 0; i < out.Rows(); i++ {
		for j := 0; j < out.Cols(); j++ {
			p := m.Row(i)
			q := b.Col(j)
			out.Set(i, j, p.Dot(&q))
		}
	}
	return out
}