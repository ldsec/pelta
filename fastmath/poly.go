package fastmath

import (
	"fmt"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Poly represents a polynomial.
type Poly struct {
	ref      *ring.Poly
	baseRing *ring.Ring
}

// NewPoly returns a zero polynomial.
func NewPoly(baseRing *ring.Ring) *Poly {
	return &Poly{baseRing.NewPoly(), baseRing}
}

func ForceInvNTT(polyNTT *PolyNTT) *Poly {
	return polyNTT.actual
}

// NewOnePoly returns a one polynomial scaled with the given factor.
func NewOnePoly(scale uint64, baseRing *ring.Ring) *Poly {
	p := NewPoly(baseRing)
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

// IsZero returns true if all the coefficients of this polynomial are zero.
func (p *Poly) IsZero() bool {
	for _, coeff := range p.ref.Coeffs[0] {
		if coeff != 0 {
			return false
		}
	}
	return true
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

// Add adds two polynomials.
func (p *PolyNTT) Add(q *PolyNTT) *PolyNTT {
	p.actual.Add(q.actual)
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

// Copy returns a copy of this polynomial.
func (p *PolyNTT) Copy() *PolyNTT {
	return &PolyNTT{&Poly{p.actual.ref.CopyNew(), p.actual.baseRing}}
}

// Coeffs returns a view into the coefficients of this polynomial.
func (p *PolyNTT) Coeffs() *IntVec {
	return &IntVec{
		size:     p.actual.N(),
		polys:    []*Poly{p.actual},
		baseRing: p.actual.baseRing,
	}
}
