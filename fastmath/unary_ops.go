package fastmath

import (
	"math"
	"math/big"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// Scale scales this polynomial with the given scalar factor.
func (p *Poly) Scale(factor uint64) *Poly {
	if p.IsUnset() {
		return p
	}
	if factor == 0 {
		return p.Zero()
	}
	p.SetDirty()
	p.baseRing.MulScalar(p.ref, factor, p.ref)
	return p
}

// Zero resets the coefficients of this polynomial to zero.
func (p *Poly) Zero() *Poly {
	p.unset = true
	p.ref.Zero()
	return p
}

// Scale scales this polynomial with the given scalar factor.
func (p *PolyNTT) Scale(factor uint64) *PolyNTT {
	p.actual.Scale(factor)
	return p
}

// Scale scales all polynomials with the given factor.
func (v *PolyVec) Scale(factor uint64) *PolyVec {
	for _, p := range v.elems {
		p.Scale(factor)
	}
	return v
}

// Scale scales all polynomials with the given factor.
func (v *PolyNTTVec) Scale(factor uint64) *PolyNTTVec {
	for _, p := range v.elems {
		p.Scale(factor)
	}
	return v
}

// Scale scales the vector by the given amount.
func (v *IntVec) Scale(factor uint64) *IntVec {
	for _, p := range v.polys {
		p.Scale(factor)
	}
	return v
}

// Scale scales the matrix by the given amount.
func (m *IntMatrix) Scale(factor uint64) *IntMatrix {
	for _, row := range m.rows {
		row.Scale(factor)
	}
	return m
}

// Neg negates this polynomial.
func (p *Poly) Neg() *Poly {
	p.baseRing.Neg(p.ref, p.ref)
	return p
}

// Neg negates this polynomial.
func (p *PolyNTT) Neg() *PolyNTT {
	p.actual.Neg()
	return p
}

// Neg negates the polynomial.
func (v *IntVec) Neg() *IntVec {
	for _, p := range v.polys {
		p.Neg()
	}
	return v
}

// Neg negates this matrix.
func (m *IntMatrix) Neg() *IntMatrix {
	for _, row := range m.rows {
		row.Neg()
	}
	return m
}

// Diag returns a diagonalization of this vector.
func (v *IntVec) Diag() *IntMatrix {
	m := NewIntMatrix(v.Size(), v.Size(), v.BaseRing())
	for i := 0; i < v.Size(); i++ {
		m.Set(i, i, v.Get(i))
	}
	return m
}

// Max returns the maximum coefficient at the given level.
func (p *Poly) Max(level int) uint64 {
	if p.IsUnset() {
		return 0
	}
	max := p.ref.Coeffs[level][0]
	for _, coeff := range p.ref.Coeffs[level][1:] {
		if coeff > max {
			max = coeff
		}
	}
	return max
}

// Max returns the maximum coefficient at the given level.
func (p *PolyNTT) Max(level int) uint64 {
	return p.actual.Max(level)
}

// Max returns the largest element of the vector.
func (v *IntVec) Max() uint64 {
	max := v.polys[0].Max(0)
	for _, p := range v.polys[1:] {
		pMax := p.Max(0)
		if pMax > max {
			max = pMax
		}
	}
	return max
}

// Max returns the largest element of the matrix.
func (v *IntMatrix) Max() uint64 {
	max := v.rows[0].Max()
	for _, r := range v.rows {
		c := r.Max()
		if c > max {
			max = c
		}
	}
	return max
}

// Max returns the largest coefficient among all the polynomials in this vector.
func (v *PolyVec) Max() uint64 {
	max := v.elems[0].Max(0)
	for _, p := range v.elems[1:] {
		c := p.Max(0)
		if c > max {
			max = c
		}
	}
	return max
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

// SumCoeffsFast returns the sum of the coefficients of this polynomial.
func (p *Poly) SumCoeffsFast(level int) uint64 {
	logN := int(math.Log2(float64(p.baseRing.N)))
	tmp := p.Copy()
	tmp2 := NewPoly(p.baseRing)
	for i := 0; i < logN; i++ {
		p.baseRing.Shift(tmp.ref, 1<<i, tmp2.ref)
		p.baseRing.AddLvl(level, tmp.ref, tmp2.ref, tmp.ref)
	}
	return tmp.ref.Coeffs[level][0]
}

// SumCoeffs returns the sum of the coefficients of this polynomial.
func (p *Poly) SumCoeffs(level int, mod *big.Int) uint64 {
	if p.IsUnset() {
		return 0
	}
	out := uint64(0)
	for _, c := range p.ref.Coeffs[level] {
		out = (out + c) % mod.Uint64()
	}
	return out
}

// Sum sums the elements of this vector.
func (v *PolyVec) Sum() *Poly {
	out := NewPoly(v.baseRing)
	for _, el := range v.elems {
		out.Add(el)
	}
	return out
}

// Sum sums the elements of this vector.
func (v *PolyNTTVec) Sum() *PolyNTT {
	out := NewPoly(v.baseRing).NTT()
	for _, el := range v.elems {
		out.Add(el)
	}
	return out
}

func (m *PolyMatrix) Sum() *Poly {
	out := NewPoly(m.baseRing)
	for _, row := range m.rows {
		rowSum := row.Sum()
		out.Add(rowSum)
	}
	return out
}

func (m *PolyNTTMatrix) Sum() *PolyNTT {
	out := NewPoly(m.baseRing).NTT()
	for _, row := range m.rows {
		rowSum := row.Sum()
		out.Add(rowSum)
	}
	return out
}

// NTT converts this polynomial to its NTT domain.
func (p *Poly) NTT() *PolyNTT {
	p.baseRing.NTT(p.ref, p.ref)
	return ForceNTT(p)
}

// NTT converts the elements of this vector to NTT space.
func (v *PolyVec) NTT() *PolyNTTVec {
	nttPolys := make([]*PolyNTT, 0, v.Size())
	for _, p := range v.elems {
		nttPolys = append(nttPolys, p.NTT())
	}
	return &PolyNTTVec{nttPolys, v.baseRing}
}

// NTT converts the elements of this matrix to NTT space.
func (m *PolyMatrix) NTT() *PolyNTTMatrix {
	nttRows := make([]*PolyNTTVec, 0, m.Rows())
	for _, r := range m.rows {
		nttRows = append(nttRows, r.NTT())
	}
	return &PolyNTTMatrix{nttRows, m.baseRing}
}

// InvNTT converts this polynomial back into its poly space.
func (p *PolyNTT) InvNTT() *Poly {
	p.actual.baseRing.InvNTT(p.actual.ref, p.actual.ref)
	return p.actual
}

// InvNTT converts the polynomials of this vector back into their poly space.
func (v *PolyNTTVec) InvNTT() *PolyVec {
	polys := make([]*Poly, 0, v.Size())
	for _, p := range v.elems {
		polys = append(polys, p.InvNTT())
	}
	return &PolyVec{polys, v.baseRing}
}

// InvNTT converts the polynomials of this matrix back into their poly space.
func (m *PolyNTTMatrix) InvNTT() *PolyMatrix {
	polyRows := make([]*PolyVec, 0, m.Rows())
	for _, r := range m.rows {
		polyRows = append(polyRows, r.InvNTT())
	}
	return &PolyMatrix{polyRows, m.baseRing}
}

// PowCoeffs takes the `exp`-th power of the coefficients modulo `mod`.
func (p *Poly) PowCoeffs(exp uint64, mod uint64) *Poly {
	if p.IsUnset() {
		return p
	}
	for i := 0; i < p.baseRing.N; i++ {
		newCoeff := ring.ModExp(p.Get(i, 0), exp, mod)
		p.Set(i, newCoeff)
	}
	return p
}

// Pow takes the exp-th power of this polynomial.
func (p *PolyNTT) Pow(exp, mod uint64) *PolyNTT {
	p.actual.PowCoeffs(exp, mod)
	return p
}

// Transposed returns the transposed version of this matrix.
func (m *IntMatrix) Transposed() *IntMatrix {
	return logging.LogShortExecution("IntMatrix.Transposed", m.SizeString(), func() interface{} {
		At := make([]uint64, m.Cols()*m.Rows())
		for row := 0; row < m.Rows(); row++ {
			for col := 0; col < m.Cols(); col++ {
				index := col*m.Rows() + row
				At[index] = m.Get(row, col)
			}
		}
		mt := NewIntMatrix(m.Cols(), m.Rows(), m.BaseRing())
		for row := 0; row < mt.Rows(); row++ {
			for col := 0; col < mt.Cols(); col++ {
				mt.Set(row, col, At[row*mt.Cols()+col])
			}
		}
		return mt
	}).(*IntMatrix)
}
