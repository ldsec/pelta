package fastmath

import "github.com/tuneinsight/lattigo/v4/ring"

type PolyMatrix struct {
	rows     []PolyVec
	baseRing *ring.Ring
}

func NewPolyMatrix(numRows, numCols int, baseRing *ring.Ring) *PolyMatrix {
	rows := make([]PolyVec, numRows)
	for i := 0; i < numRows; i++ {
		rows[i] = *NewPolyVec(numCols, baseRing)
	}
	return &PolyMatrix{rows, baseRing}
}

func (m *PolyMatrix) Rows() int {
	return len(m.rows)
}

func (m *PolyMatrix) Cols() int {
	return m.rows[0].Size()
}

func (m *PolyMatrix) Row(i int) *PolyVec {
	return &m.rows[i]
}

func (m *PolyMatrix) Get(i, j int) *Poly {
	return m.rows[i].Get(j)
}

func (m *PolyMatrix) Populate(f func(int, int) *Poly) {
	for i := 0; i < len(m.rows); i++ {
		m.rows[i].Populate(func(j int) *Poly {
			return f(i, j)
		})
	}
}

func (m *PolyMatrix) Update(f func(int, int, *Poly) *Poly) {
	for i := range m.rows {
		m.rows[i].Update(func(j int, old *Poly) *Poly {
			return f(i, j, old)
		})
	}
}

func (m *PolyMatrix) PopulateRows(f func(int) PolyVec) {
	for i := range m.rows {
		m.rows[i] = f(i)
	}
}

func (m *PolyMatrix) UpdateRows(f func(int, PolyVec) PolyVec) {
	for i, old := range m.rows {
		m.rows[i] = f(i, old)
	}
}

func (m *PolyMatrix) Sum() *Poly {
	out := NewZeroPoly(m.baseRing)
	for _, row := range m.rows {
		rowSum := row.Sum()
		out.Add(rowSum)
	}
	return out
}

func (m *PolyMatrix) Copy() *PolyMatrix {
	newRows := make([]PolyVec, 0, len(m.rows))
	for _, row := range m.rows {
		newRows = append(newRows, *row.Copy())
	}
	return &PolyMatrix{newRows, m.baseRing}
}

func (m *PolyMatrix) AllRows(pred func(int, *PolyVec) bool) bool {
	for i, row := range m.rows {
		if !pred(i, &row) {
			return false
		}
	}
	return true
}

func (m *PolyMatrix) NTT() *PolyNTTMatrix {
	nttRows := make([]PolyNTTVec, 0, m.Rows())
	for _, r := range m.rows {
		nttRows = append(nttRows, *r.NTT())
	}
	return &PolyNTTMatrix{nttRows, m.baseRing}
}

func (m *PolyMatrix) Eq(b *PolyMatrix) bool {
	if m.Rows() != b.Rows() || m.Cols() != b.Cols() {
		return false
	}
	for i, r := range m.rows {
		if !r.Eq(b.Row(i)) {
			return false
		}
	}
	return true
}

type PolyNTTMatrix struct {
	rows     []PolyNTTVec
	baseRing *ring.Ring
}

func (m *PolyNTTMatrix) Rows() int {
	return len(m.rows)
}

func (m *PolyNTTMatrix) Cols() int {
	return m.rows[0].Size()
}

func (m *PolyNTTMatrix) Row(i int) *PolyNTTVec {
	return &m.rows[i]
}

func (m *PolyNTTMatrix) ColCopy(i int) *PolyNTTVec {
	col := NewPolyVec(m.Rows(), m.baseRing).NTT()
	for j := 0; j < m.Rows(); j++ {
		col.Set(j, m.Get(i, j))
	}
	return col
}

func (m *PolyNTTMatrix) Get(i, j int) *PolyNTT {
	return m.rows[i].Get(j)
}

func (m *PolyNTTMatrix) Set(i, j int, val *PolyNTT) {
	m.rows[i].Set(j, val)
}

func (m *PolyNTTMatrix) Populate(f func(int, int) *PolyNTT) {
	for i := 0; i < len(m.rows); i++ {
		m.rows[i].Populate(func(j int) *PolyNTT {
			return f(i, j)
		})
	}
}

func (m *PolyNTTMatrix) Update(f func(int, int, *PolyNTT) *PolyNTT) {
	for i := range m.rows {
		m.rows[i].Update(func(j int, old *PolyNTT) *PolyNTT {
			return f(i, j, old)
		})
	}
}

func (m *PolyNTTMatrix) PopulateRows(f func(int) PolyNTTVec) {
	for i := range m.rows {
		m.rows[i] = f(i)
	}
}

func (m *PolyNTTMatrix) UpdateRows(f func(int, PolyNTTVec) PolyNTTVec) {
	for i, old := range m.rows {
		m.rows[i] = f(i, old)
	}
}

func (m *PolyNTTMatrix) Sum() *PolyNTT {
	out := NewZeroPoly(m.baseRing).NTT()
	for _, row := range m.rows {
		rowSum := row.Sum()
		out.Add(rowSum)
	}
	return out
}

func (m *PolyNTTMatrix) MulVec(b *PolyNTTVec) *PolyNTTVec {
	out := NewPolyVec(m.Rows(), m.baseRing).NTT()
	for i, row := range m.rows {
		prod := row.Dot(b)
		out.Set(i, prod)
	}
	return out
}

// MulMat performs a matrix-matrix multiplication.
func (m *PolyNTTMatrix) MulMat(b *PolyNTTMatrix) *PolyNTTMatrix {
	if m.Cols() != b.Rows() {
		panic("IntMatrix.MulMat sizes incorrect")
	}
	out := NewPolyMatrix(m.Rows(), b.Cols(), m.baseRing).NTT()
	for i := 0; i < out.Rows(); i++ {
		for j := 0; j < out.Cols(); j++ {
			p := m.Row(i)
			q := b.ColCopy(j)
			out.Set(i, j, p.Dot(q))
		}
	}
	return out
}

func (m *PolyNTTMatrix) Copy() *PolyNTTMatrix {
	newRows := make([]PolyNTTVec, 0, len(m.rows))
	for _, row := range m.rows {
		newRows = append(newRows, *row.Copy())
	}
	return &PolyNTTMatrix{newRows, m.baseRing}
}

func (m *PolyNTTMatrix) AllRows(pred func(int, *PolyNTTVec) bool) bool {
	for i, row := range m.rows {
		if !pred(i, &row) {
			return false
		}
	}
	return true
}

func (m *PolyNTTMatrix) All(pred func(int, int, *PolyNTT) bool) bool {
	for i, row := range m.rows {
		if !row.All(func(j int, el *PolyNTT) bool {
			return pred(i, j, el)
		}) {
			return false
		}
	}
	return true
}

func (m *PolyNTTMatrix) ExtendCols(a *PolyNTTMatrix) {
	for i, r := range m.rows {
		for _, el := range a.Row(i).elems {
			r.Append(el)
		}
	}
}

func (m *PolyNTTMatrix) ExtendRows(a *PolyNTTMatrix) {
	for _, r := range a.rows {
		m.rows = append(m.rows, r)
	}
}

func (m *PolyNTTMatrix) InvNTT() *PolyMatrix {
	polyRows := make([]PolyVec, 0, m.Rows())
	for _, r := range m.rows {
		polyRows = append(polyRows, *r.InvNTT())
	}
	return &PolyMatrix{polyRows, m.baseRing}
}

func (m *PolyNTTMatrix) Eq(b *PolyNTTMatrix) bool {
	if m.Rows() != b.Rows() || m.Cols() != b.Cols() {
		return false
	}
	for i, r := range m.rows {
		if !r.Eq(b.Row(i)) {
			return false
		}
	}
	return true
}

// ToIntMatrix converts a polynomial matrix of dimensions m x n to an integer matrix of dimensions m x nd.
func (m *PolyNTTMatrix) ToIntMatrix() *IntMatrix {
	d := m.baseRing.N
	cols := m.Cols() * d
	im := NewIntMatrix(m.Rows(), cols, m.baseRing)
	im.PopulateRows(func(i int) *IntVec {
		return m.Row(i).ToIntVec()
	})
	return im
}
