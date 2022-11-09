package fastmath

import (
	"fmt"
	"math"
	"strings"

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

// N returns the degree of the polynomial.
func (p *Poly) N() int {
	return p.ref.N()
}

// Coeffs returns a view into the coefficients of this polynomial.
func (p *Poly) Coeffs() IntVec {
	return IntVec{size: p.N(), polys: []Poly{*p}, baseRing: p.baseRing}
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

// NTT converts this polynomial to its NTT domain.
func (p *Poly) NTT() *Poly {
	p.baseRing.NTT(p.ref, p.ref)
	return p
}

// InvNTT converts this polynomial back to its poly domain.
func (p *Poly) InvNTT() *Poly {
	p.baseRing.InvNTT(p.ref, p.ref)
	return p
}

// Neg negates this polynomial.
func (p *Poly) Neg() *Poly {
	p.baseRing.Neg(p.ref, p.ref)
	return p
}

// Scale scales this polynomial with the given scalar factor.
func (p *Poly) Scale(factor uint64) *Poly {
	p.baseRing.MulScalar(p.ref, factor, p.ref)
	return p
}

// Add adds two polynomials.
func (p *Poly) Add(q *Poly) *Poly {
	p.baseRing.Add(p.ref, q.ref, p.ref)
	return p
}

// MulCoeffs multiplies the coefficients of two polynomials.
func (p *Poly) MulCoeffs(q *Poly) *Poly {
	p.baseRing.MulCoeffs(p.ref, q.ref, p.ref)
	return p
}

// Pow takes the `exp`-th power of the coefficients modulo `mod`.
func (p *Poly) PowModCoeffs(exp uint64, mod uint64) *Poly {
	for i := 0; i < p.baseRing.N; i++ {
		newCoeff := ring.ModExp(p.Get(i, 0), exp, mod)
		p.Set(i, newCoeff)
	}
	return p
}

// Reduce performs the appropriate coefficient reductions over all the levels.
func (p *Poly) Reduce() *Poly {
	p.baseRing.Reduce(p.ref, p.ref)
	return p
}

func (p *Poly) Eq(q *Poly) bool {
	return p.baseRing.Equal(p.ref, q.ref)
}

func (p *Poly) EqLevel(level int, q *Poly) bool {
	return p.baseRing.EqualLvl(level, p.ref, q.ref)
}

func (p *Poly) String() string {
	s := "Poly{"
	coeffStrings := make([]string, 0)
	for _, c := range p.ref.Coeffs[0][:10] {
		coeffStrings = append(coeffStrings, fmt.Sprintf("%d", c))
	}
	return s + strings.Join(coeffStrings, ",") + "}"
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

func (v *PolyVec) Populate(f func(int) Poly) {
	for i := 0; i < len(v.elems); i++ {
		v.elems[i] = f(i)
	}
}

func (v *PolyVec) Sum() Poly {
	out := NewZeroPoly(v.baseRing)
	for _, el := range v.elems {
		out.Add(&el)
	}
	return out
}

type PolyMatrix struct {
	rows     []PolyVec
	baseRing *ring.Ring
}

func NewPolyMatrix(numRows, numCols int, baseRing *ring.Ring) PolyMatrix {
	rows := make([]PolyVec, numRows)
	for i := 0; i < numRows; i++ {
		rows[i] = NewPolyVec(numCols, baseRing)
	}
	return PolyMatrix{rows, baseRing}
}

func (m *PolyMatrix) Populate(f func(int, int) Poly) {
	for i := 0; i < len(m.rows); i++ {
		m.rows[i].Populate(func(j int) Poly {
			return f(i, j)
		})
	}
}

func (m *PolyMatrix) PopulateRows(f func(int) PolyVec) {
	for i := 0; i < len(m.rows); i++ {
		m.rows[i] = f(i)
	}
}