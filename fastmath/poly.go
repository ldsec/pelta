package fastmath

import (
	"fmt"
	"math"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Poly represents a polynomial.
type Poly struct {
	ref      *ring.Poly
	baseRing *ring.Ring
}

// NewZeroPoly returns a zero polynomial.
func NewZeroPoly(baseRing *ring.Ring) *Poly {
	return &Poly{baseRing.NewPoly(), baseRing}
}

// NewOnePoly returns a one polynomial scaled with the given factor.
func NewOnePoly(scale uint64, baseRing *ring.Ring) *Poly {
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
func (p *Poly) Coeffs() *IntVec {
	return &IntVec{size: p.N(), polys: []*Poly{p}, baseRing: p.baseRing}
}

// SumCoeffs returns the sum of the coefficients of this polynomial.
func (p *Poly) SumCoeffs(level int) uint64 {
	logN := int(math.Log2(float64(p.baseRing.N)))
	tmp := p.Copy()
	tmp2 := NewZeroPoly(p.baseRing)
	for i := 0; i < logN; i++ {
		p.baseRing.Shift(tmp.ref, 1<<i, tmp2.ref)
		p.baseRing.Add(tmp.ref, tmp2.ref, tmp.ref)
	}
	return tmp.ref.Coeffs[level][0]
}

// SumCoeffsLimited returns the sum of the first `limit` coefficients of this polynomial.
func (p *Poly) SumCoeffsLimited(level, limit int, mod *big.Int) uint64 {
	if limit >= p.N() {
		return p.SumCoeffs(level)
	}
	out := big.NewInt(0)
	for i := 0; i < limit; i++ {
		out.Add(out, big.NewInt(int64(p.Get(i, level)))).Mod(out, mod)
	}
	return out.Uint64()
}

// MaxCoeff returns the maximum coefficient at the given level.
func (p *Poly) MaxCoeff(level int) uint64 {
	max := p.ref.Coeffs[level][0]
	for _, coeff := range p.ref.Coeffs[level][1:] {
		if coeff > max {
			max = coeff
		}
	}
	return max
}

// IsZero returns true if all the coefficients of this polynomial are zero.
func (p *Poly) IsZero() bool {
	for _, coeff := range p.ref.Coeffs[0] {
		if coeff != 0 {
			return false
		}
	}
	return true
}

// NTT converts this polynomial to its NTT domain.
func (p *Poly) NTT() *PolyNTT {
	c := p.Copy()
	p.baseRing.NTT(c.ref, c.ref)
	return &PolyNTT{c}
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

// PowCoeffs takes the `exp`-th power of the coefficients modulo `mod`.
func (p *Poly) PowCoeffs(exp uint64, mod uint64) *Poly {
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

// Eq checks the equality of the coefficients of the two rings.
func (p *Poly) Eq(q *Poly) bool {
	return p.baseRing.Equal(p.ref, q.ref)
}

// EqLevel checks the equality of the coefficients of the two rings at the given level.
func (p *Poly) EqLevel(level int, q *Poly) bool {
	return p.baseRing.EqualLvl(level, p.ref, q.ref)
}

// String returns a string representation of the polynomial.
func (p *Poly) String() string {
	s := fmt.Sprintf("Poly{")
	coeffStrings := make([]string, 0)
	for _, c := range p.ref.Coeffs[0][:10] {
		coeffStrings = append(coeffStrings, fmt.Sprintf("%d", c))
	}
	return s + strings.Join(coeffStrings, ",") + "}"
}

// Copy returns a copy of this polynomial.
func (p *Poly) Copy() *Poly {
	return &Poly{p.ref.CopyNew(), p.baseRing}
}

// PolyNTT represents a polynomial in the NTT domain.
type PolyNTT struct {
	actual *Poly
}

func ForceNTT(p *Poly) *PolyNTT {
	return &PolyNTT{p}
}

func (p *PolyNTT) Get(i, level int) uint64 {
	return p.actual.Get(i, level)
}

// Mul multiplies two polynomials.
func (p *PolyNTT) Mul(q *PolyNTT) *PolyNTT {
	p.actual.MulCoeffs(q.actual)
	return p
}

// InvNTT converts this polynomial back into its poly space.
func (p *PolyNTT) InvNTT() *Poly {
	p.actual.baseRing.InvNTT(p.actual.ref, p.actual.ref)
	return p.actual
}

// Pow takes the exp-th power of this polynomial.
func (p *PolyNTT) Pow(exp, mod uint64) *PolyNTT {
	p.actual.PowCoeffs(exp, mod)
	return p
}

// Add adds two polynomials.
func (p *PolyNTT) Add(q *PolyNTT) *PolyNTT {
	p.actual.Add(q.actual)
	return p
}

// Scale scales this polynomial with the given scalar factor.
func (p *PolyNTT) Scale(factor uint64) *PolyNTT {
	p.actual.Scale(factor)
	return p
}

// String returns a string representation of this polynomial.
func (p *PolyNTT) String() string {
	s := fmt.Sprintf("PolyNTT{")
	coeffStrings := make([]string, 0)
	for _, c := range p.actual.ref.Coeffs[0][:10] {
		coeffStrings = append(coeffStrings, fmt.Sprintf("%d", c))
	}
	return s + strings.Join(coeffStrings, ",") + "}"
}

// Neg negates this polynomial.
func (p *PolyNTT) Neg() *PolyNTT {
	p.actual.Neg()
	return p
}

// Copy returns a copy of this polynomial.
func (p *PolyNTT) Copy() *PolyNTT {
	return &PolyNTT{&Poly{p.actual.ref.CopyNew(), p.actual.baseRing}}
}

// Eq checks whether these two polynomials are equal.
func (p *PolyNTT) Eq(q *PolyNTT) bool {
	return p.actual.baseRing.EqualLvl(0, p.actual.ref, q.actual.ref)
}

// Coeffs returns a view into the coefficients of this polynomial.
func (p *PolyNTT) Coeffs() *IntVec {
	return &IntVec{
		size:     p.actual.N(),
		polys:    []*Poly{p.actual},
		baseRing: p.actual.baseRing,
	}
}

// MaxCoeff returns the maximum coefficient at the given level.
func (p *PolyNTT) MaxCoeff(level int) uint64 {
	max := p.actual.ref.Coeffs[level][0]
	for _, coeff := range p.actual.ref.Coeffs[level][1:] {
		if coeff > max {
			max = coeff
		}
	}
	return max
}
