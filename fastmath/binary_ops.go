package fastmath

import (
	"math/big"
)

// Add adds two polynomials.
func (p *Poly) Add(q *Poly) *Poly {
	// if q.IsZero() {
	// 	return p
	// }
	p.unset = false
	p.baseRing.Add(p.ref, q.ref, p.ref)
	return p
}

// MulCoeffs multiplies the coefficients of two polynomials.
func (p *Poly) MulCoeffs(q *Poly) *Poly {
	// if q.IsZero() {
	// 	return p.Zero()
	// }
	p.unset = false
	p.baseRing.MulCoeffs(p.ref, q.ref, p.ref)
	return p
}

// MulCoeffsAndAdd multiplies the coefficients of two polynomials and adds it to `out`.
func (p *Poly) MulCoeffsAndAdd(q *Poly, out *Poly) *Poly {
	// if q.IsZero() {
	// 	return p.Zero()
	// }
	out.unset = false
	p.baseRing.MulCoeffsAndAdd(p.ref, q.ref, out.ref)
	return p
}

// Mul multiplies two polynomials.
func (p *PolyNTT) Mul(q *PolyNTT) *PolyNTT {
	p.actual.MulCoeffs(q.actual)
	return p
}

// Add adds two polynomials.
func (p *PolyNTT) Add(q *PolyNTT) *PolyNTT {
	p.actual.Add(q.actual)
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
		a.MulCoeffsAndAdd(b, out)
	}
}

// Dot returns the dot product of the given two vectors.
func (v *IntVec) Dot(r *IntVec) uint64 {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	if v.mod == nil {
		v.mod = v.baseRing.ModulusAtLevel[0]
	}
	preSum := NewPoly(v.baseRing)
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		a.MulCoeffsAndAdd(b, preSum)
	}
	return preSum.SumCoeffs(0, v.mod)
}

func (v *PolyNTTVec) Dot(r *PolyNTTVec) *PolyNTT {
	out := NewPoly(v.baseRing)
	for i, p := range v.elems {
		a := p.actual
		b := r.Get(i).actual
		a.MulCoeffsAndAdd(b, out)
	}
	return ForceNTT(out)
}

// Hadamard performs coefficient-wise multiplication.
func (v *IntVec) Hadamard(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	for i := 0; i < len(v.polys); i++ {
		a := v.polys[i]
		b := r.polys[i]
		a.MulCoeffs(b)
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

// EqLevel checks the equality of the coefficients of the two rings up to the given level.
func (p *Poly) EqLevel(level int, q *Poly) bool {
	return p.baseRing.EqualLvl(level, p.ref, q.ref)
}

// Eq checks whether these two polynomials are equal.
func (p *PolyNTT) Eq(q *PolyNTT) bool {
	// return p.actual.baseRing.EqualLvl(0, p.actual.ref, q.actual.ref)
	return p.actual.baseRing.Equal(p.actual.ref, p.actual.ref)
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
		v.baseRing.Reduce(v.polys[i].ref, v.polys[i].ref)
		v.baseRing.Reduce(r.polys[i].ref, r.polys[i].ref)
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
		out.Set(i, dotResult.SumCoeffs(0, m.mod))
		dotResult.Zero()
	}
	return out
}

// MulVec performs a matrix-vector multiplication.
func (m *IntMatrix) MulVec2(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("IntMatrix.MulVec sizes incorrect")
	}
	out := NewIntVec(m.Rows(), m.baseRing)
	for i, row := range m.rows {
		dotResult := big.NewInt(0)
		for j := 0; j < m.Cols(); j++ {
			mult := big.NewInt(0).Mul(big.NewInt(int64(v.Get(j))), big.NewInt(int64(row.Get(j))))
			dotResult.Add(dotResult, mult)
			dotResult.Mod(dotResult, m.mod)
		}
		out.Set(i, dotResult.Mod(dotResult, m.mod).Uint64())
	}
	return out
}

// MulVecTransposed performs a matrix-vector multiplication with the transpose of this matrix.
func (m *IntMatrix) MulVecTransposed(v *IntVec) *IntVec {
	if m.Rows() != v.Size() {
		panic("IntMatrix.MulVecTransposed sizes incorrect")
	}
	out := NewIntVec(m.Cols(), m.baseRing)
	for j := 0; j < m.Cols(); j++ {
		dotResult := big.NewInt(0)
		for i := 0; i < m.Rows(); i++ {
			mult := big.NewInt(0).Mul(big.NewInt(int64(v.Get(i))), big.NewInt(int64(m.Get(i, j))))
			dotResult.Add(dotResult, mult)
			dotResult.Mod(dotResult, m.mod)
		}
		out.Set(j, dotResult.Mod(dotResult, m.mod).Uint64())
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
	modUint := m.mod.Uint64()
	for i := 0; i < out.Rows(); i++ {
		for j := 0; j < out.Cols(); j++ {
			p := m.RowView(i)
			dotResult := uint64(0)
			for k := 0; k < p.Size(); k++ {
				dotResult = (dotResult + p.Get(k)*b.Get(k, j)) % modUint
			}
			out.Set(i, j, dotResult)
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
