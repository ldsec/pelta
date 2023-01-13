package fastmath

import (
	"fmt"
	"math"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Poly represents a polynomial.
type Poly struct {
	ref      *ring.Poly
	baseRing *ring.Ring
	unset    bool
}

// NewPoly returns a zero polynomial.
func NewPoly(baseRing *ring.Ring) *Poly {
	return &Poly{baseRing.NewPoly(), baseRing, true}
}

func ForceInvNTT(polyNTT *PolyNTT) *Poly {
	return polyNTT.actual
}

// NewOnePoly returns a one polynomial scaled with the given factor.
func NewOnePoly(scale uint64, baseRing *ring.Ring) *Poly {
	p := NewPoly(baseRing)
	p.SetForce(0, scale)
	return p
}

// NewOnePoly returns a one polynomial scaled with the given factor.
func NewOnePolyLevels(scales []uint64, baseRing *ring.Ring) *Poly {
	p := NewPoly(baseRing)
	p.SetCoeff(0, scales)
	return p
}

// SetCoeff sets the coefficients of this polynomial at every level to the given values.
func (p *Poly) SetCoeff(index int, coeff Coeff) {
	if !coeff.IsZero() {
		p.SetDirty()
	}
	for level := 0; level < len(p.ref.Coeffs); level++ {
		p.SetLevel(index, level, coeff[level])
	}
}

// SetForce sets the coefficient of this polynomial at every level to the given value.
func (p *Poly) SetForce(index int, value uint64) {
	if value != 0 {
		p.SetDirty()
	}
	for level := 0; level < len(p.ref.Coeffs); level++ {
		p.SetLevel(index, level, value)
	}
}

// SetLevel sets the coefficient of this polynomial at the given level to the given value.
func (p *Poly) SetLevel(index, level int, value uint64) {
	if value != 0 {
		p.SetDirty()
	}
	p.ref.Coeffs[level][index] = value
}

// GetCoeff returns the coefficient of this polynomial.
func (p *Poly) GetCoeff(index int) Coeff {
	coeffs := make([]uint64, len(p.ref.Coeffs))
	for lvl := 0; lvl < len(p.ref.Coeffs); lvl++ {
		coeffs[lvl] = p.ref.Coeffs[lvl][index]
	}
	return coeffs
}

// GetLevel returns the coefficient of this polynomial at the given level.
func (p *Poly) GetLevel(index, level int) uint64 {
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

func (p *Poly) IsUnset() bool {
	return p.unset
}

// IsZero returns true if all the coefficients of this polynomial are zero.
func (p *Poly) IsZero() bool {
	if p.IsUnset() {
		return true
	}
	for _, coeff := range p.ref.Coeffs[0] {
		if coeff != 0 {
			return false
		}
	}
	return true
}

func (p *Poly) SetDirty() {
	p.unset = false
}

// Reduce performs the appropriate coefficient reductions over all the levels.
func (p *Poly) Reduce() *Poly {
	p.baseRing.Reduce(p.ref, p.ref)
	return p
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
	return &Poly{p.ref.CopyNew(), p.baseRing, p.unset}
}

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

// Scale scales this polynomial with the given coefficient.
func (p *Poly) ScaleCoeff(factors Coeff) *Poly {
	if p.IsUnset() {
		return p
	}
	if factors.IsZero() {
		return p.Zero()
	}
	p.SetDirty()
	p.MulCoeffs(factors.ExtendAsPoly(p.baseRing))
	return p
}

// Zero resets the coefficients of this polynomial to zero.
func (p *Poly) Zero() *Poly {
	p.unset = true
	p.ref.Zero()
	return p
}

// Neg negates this polynomial.
func (p *Poly) Neg() *Poly {
	p.baseRing.Neg(p.ref, p.ref)
	return p
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

// SumCoeffs returns the sum of the coefficients of this polynomial.
func (p *Poly) SumCoeffs() Coeff {
	logN := int(math.Log2(float64(p.baseRing.N)))
	tmp := p.Copy()
	tmp2 := NewPoly(p.baseRing)
	for i := 0; i < logN; i++ {
		p.baseRing.Shift(tmp.ref, 1<<i, tmp2.ref)
		p.baseRing.Add(tmp.ref, tmp2.ref, tmp.ref)
	}
	return tmp.GetCoeff(0)
}

// NTT converts this polynomial to its NTT domain.
func (p *Poly) NTT() *PolyNTT {
	p.baseRing.NTT(p.ref, p.ref)
	return ForceNTT(p)
}

// PowCoeffs takes the `exp`-th power of the coefficients.
func (p *Poly) PowCoeffs(exp uint) *Poly {
	if p.IsUnset() {
		return p
	}
	mul := p.Copy()
	for i := 0; i < int(exp); i++ {
		p.MulCoeffs(mul)
	}
	return p
}

// Add adds two polynomials.
func (p *Poly) Add(q *Poly) *Poly {
	if q.IsUnset() {
		return p
	}
	p.SetDirty()
	p.baseRing.Add(p.ref, q.ref, p.ref)
	return p
}

// MulCoeffs multiplies the coefficients of two polynomials.
func (p *Poly) MulCoeffs(q *Poly) *Poly {
	if q.IsUnset() || p.IsUnset() {
		return p.Zero()
	}
	p.SetDirty()
	p.baseRing.MulCoeffs(p.ref, q.ref, p.ref)
	return p
}

// MulCoeffsAndAdd multiplies the coefficients of two polynomials and adds it to `out`.
func (p *Poly) MulCoeffsAndAdd(q *Poly, out *Poly) {
	if q.IsUnset() || p.IsUnset() {
		return
	}
	out.SetDirty()
	p.baseRing.MulCoeffsAndAdd(p.ref, q.ref, out.ref)
}

// Eq checks the equality of the coefficients of the two rings across all levels.
func (p *Poly) Eq(q *Poly) bool {
	return p.baseRing.Equal(p.ref, q.ref)
}

// EqLevel checks the equality of the coefficients of the two rings up to the given level.
func (p *Poly) EqLevel(level int, q *Poly) bool {
	return p.baseRing.EqualLvl(level, p.ref, q.ref)
}

// PolyNTT represents a polynomial in the NTT domain.
type PolyNTT struct {
	actual *Poly
}

func ForceNTT(p *Poly) *PolyNTT {
	return &PolyNTT{p}
}

func (p *PolyNTT) GetLevel(i, level int) uint64 {
	return p.actual.GetLevel(i, level)
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

// Copy returns a copy of this polynomial.
func (p *PolyNTT) Copy() *PolyNTT {
	return &PolyNTT{&Poly{p.actual.ref.CopyNew(), p.actual.baseRing, p.actual.unset}}
}

// Coeffs returns a view into the coefficients of this polynomial.
func (p *PolyNTT) Coeffs() *IntVec {
	return &IntVec{
		size:     p.actual.N(),
		polys:    []*Poly{p.actual},
		baseRing: p.actual.baseRing,
	}
}

// Pow takes the exp-th power of this polynomial.
func (p *PolyNTT) Pow(exp uint) *PolyNTT {
	p.actual.PowCoeffs(exp)
	return p
}

// Scale scales this polynomial with the given scalar factor.
func (p *PolyNTT) Scale(factor uint64) *PolyNTT {
	p.actual.Scale(factor)
	return p
}

// Neg negates this polynomial.
func (p *PolyNTT) Neg() *PolyNTT {
	p.actual.Neg()
	return p
}

// Max returns the maximum coefficient at the given level.
func (p *PolyNTT) Max(level int) uint64 {
	return p.actual.Max(level)
}

// InvNTT converts this polynomial back into its poly space.
func (p *PolyNTT) InvNTT() *Poly {
	p.actual.baseRing.InvNTT(p.actual.ref, p.actual.ref)
	return p.actual
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

// Eq checks whether these two polynomials are equal.
func (p *PolyNTT) Eq(q *PolyNTT) bool {
	// return p.actual.baseRing.EqualLvl(0, p.actual.ref, q.actual.ref)
	return p.actual.baseRing.Equal(p.actual.ref, p.actual.ref)
}
