package fastmath

// Add adds two polynomials.
func (p *Poly) Add(q *Poly) *Poly {
	p.baseRing.Add(p.ref, q.ref, p.ref)
	return p
}

// Add adds two integer vectors.
func (v *IntVec) Add(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Add sizes do not match")
	}
	for i, p := range v.polys {
		p.Add(r.polys[i])
	}
	return v
}

func (v *PolyVec) Add(b *PolyVec) *PolyVec {
	for i, p := range v.elems {
		p.Add(b.Get(i))
	}
	return v
}

func (v *PolyNTTVec) Add(b *PolyNTTVec) *PolyNTTVec {
	for i, p := range v.elems {
		p.Add(b.Get(i))
	}
	return v
}

func (v *PolyVec) AddAll(q *Poly) *PolyVec {
	for _, p := range v.elems {
		p.Add(q)
	}
	return v
}

func (v *PolyNTTVec) AddAll(q *PolyNTT) *PolyNTTVec {
	for _, p := range v.elems {
		p.Add(q)
	}
	return v
}

// MulCoeffs multiplies the coefficients of two polynomials.
func (p *Poly) MulCoeffs(q *Poly) *Poly {
	p.baseRing.MulCoeffs(p.ref, q.ref, p.ref)
	return p
}
func (v *PolyNTTVec) Mul(r *PolyNTTVec) *PolyNTTVec {
	for i, p := range v.elems {
		p.Mul(r.Get(i))
	}
	return v
}

func (v *PolyNTTVec) MulAll(r *PolyNTT) *PolyNTTVec {
	for _, p := range v.elems {
		p.Mul(r)
	}
	return v
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

// Dot returns the dot product of the given two vectors.
func (v *IntVec) Dot(r *IntVec) uint64 {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	preSum := NewPoly(v.baseRing)
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		v.baseRing.MulCoeffsAndAdd(a.ref, b.ref, preSum.ref)
	}
	return preSum.SumCoeffsLimited(0, v.size, v.mod)
}

func (v *PolyNTTVec) Dot(r *PolyNTTVec) *PolyNTT {
	out := NewPoly(v.baseRing).NTT()
	for i, p := range v.elems {
		a := p.actual.ref
		b := r.Get(i).actual.ref
		v.baseRing.MulCoeffsAndAdd(a, b, out.actual.ref)
	}
	return out
}

// Hadamard performs coefficient-wise multiplication.
func (v *IntVec) Hadamard(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		v.baseRing.MulCoeffs(a.ref, b.ref, a.ref)
	}
	return v
}

// Hadamard performs coefficient-wise multiplication.
func (m *IntMatrix) Hadamard(b *IntMatrix) *IntMatrix {
	for i, r := range m.rows {
		r.Hadamard(b.RowView(i))
	}
	return m
}

// Eq checks the equality of the coefficients of the two rings.
func (p *Poly) Eq(q *Poly) bool {
	return p.baseRing.Equal(p.ref, q.ref)
}

// EqLevel checks the equality of the coefficients of the two rings at the given level.
func (p *Poly) EqLevel(level int, q *Poly) bool {
	return p.baseRing.EqualLvl(level, p.ref, q.ref)
}

// Eq checks whether these two polynomials are equal.
func (p *PolyNTT) Eq(q *PolyNTT) bool {
	return p.actual.baseRing.EqualLvl(0, p.actual.ref, q.actual.ref)
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
		if !v.polys[i].EqLevel(0, r.polys[i]) {
			return false
		}
	}
	return true
}

func (v *PolyNTTVec) Eq(r *PolyNTTVec) bool {
	for i, p := range v.elems {
		if !p.Eq(r.Get(i)) {
			return false
		}
	}
	return true
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

// MulVec performs a matrix-vector multiplication.
func (m *IntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	out := NewIntVec(m.Rows(), m.baseRing)
	dotResult := NewPoly(m.baseRing)
	for i, row := range m.rows {
		row.MulAddElems(v, dotResult)
		out.Set(i, dotResult.SumCoeffsLimited(0, row.Size(), m.mod))
		dotResult.ref.Zero()
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
