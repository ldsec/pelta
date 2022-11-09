package fastmath

import "github.com/tuneinsight/lattigo/v4/ring"

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