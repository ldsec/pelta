package math

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
	"strings"
)

type Polynomial struct {
	Ref      *ring.Poly
	BaseRing *ring.Ring
}

func NewZeroPolynomial(baseRing *ring.Ring) Polynomial {
	return Polynomial{baseRing.NewPoly(), baseRing}
}

func NewOnePolynomial(baseRing *ring.Ring) Polynomial {
	poly := Polynomial{baseRing.NewPoly(), baseRing}
	return poly.One().(Polynomial)
}

// Copy copies the polynomial.
func (p Polynomial) Copy() RingElement {
	return Polynomial{p.Ref.CopyNew(), p.BaseRing}
}

// Add adds a polynomial to this polynomial, and returns itself.
func (p Polynomial) Add(q RingElement) RingElement {
	p.InvNTT()
	q.(Polynomial).InvNTT()
	p.BaseRing.Add(p.Ref, q.(Polynomial).Ref, p.Ref)
	//p.BaseRing.Reduce(p.Ref, p.Ref)
	return p
}

// Sub subtracts a polynomial from this polynomial, and returns itself.
func (p Polynomial) Sub(q RingElement) RingElement {
	return p.Add(q.Copy().Neg())
}

// Mul multiplies a polynomial with this polynomial, and returns itself.
func (p Polynomial) Mul(q RingElement) RingElement {
	p.NTT()
	q.(Polynomial).NTT()
	p.BaseRing.MulCoeffs(p.Ref, q.(Polynomial).Ref, p.Ref)
	//p.InvNTT()
	//q.(Polynomial).InvNTT()
	return p
}

func (p Polynomial) MulAdd(q RingElement, out RingElement) {
	p.NTT()
	q.(Polynomial).NTT()
	out.(Polynomial).NTT()
	p.BaseRing.MulCoeffsAndAdd(p.Ref, q.(Polynomial).Ref, out.(Polynomial).Ref)
	//p.InvNTT()
	//q.(Polynomial).InvNTT()
	//out.(Polynomial).InvNTT()
}

func (p Polynomial) Neg() RingElement {
	p.BaseRing.Neg(p.Ref, p.Ref)
	return p
}

// Scale scales the coefficients of the polynomial in-place.
// p => p = c*p
func (p Polynomial) Scale(factor uint64) RingElement {
	p.BaseRing.MulScalar(p.Ref, factor, p.Ref)
	return p
}

// Pow takes the exponentiation of this polynomial, returning itself.
// Warning: Unoptimized & slow
func (p Polynomial) Pow(exp uint64) RingElement {
	if exp == 0 {
		return NewOnePolynomial(p.BaseRing)
	}
	out := p.Copy()
	for i := uint64(1); i < exp; i++ {
		out.Mul(p)
	}
	//p.BaseRing.Reduce(out.(Polynomial).Ref, out.(Polynomial).Ref)
	p.Ref.CopyValues(out.(Polynomial).Ref)
	return p
}

func (p Polynomial) Eq(el RingElement) bool {
	p.InvNTT()
	p.BaseRing.Reduce(p.Ref, p.Ref)
	el.(Polynomial).InvNTT()
	el.(Polynomial).BaseRing.Reduce(el.(Polynomial).Ref, el.(Polynomial).Ref)
	for i, c1 := range p.Ref.Coeffs[0] {
		if c1 != el.(Polynomial).Ref.Coeffs[0][i] {
			return false
		}
	}
	return true
}

// Zero sets all the coefficients of this polynomial to zero.
func (p Polynomial) Zero() RingElement {
	p.Ref.Zero()
	return p
}

// One sets the first coefficient of this polynomial to one and the remaining to zero.
func (p Polynomial) One() RingElement {
	p.Zero()
	p.SetCoefficient(0, 1)
	return p
}

// NTT converts into the NTT space in-place.
// p => p = NTT(p)
func (p Polynomial) NTT() Polynomial {
	if !p.Ref.IsNTT {
		p.Ref.IsNTT = true
		p.BaseRing.NTT(p.Ref, p.Ref)
	}
	return p
}

// InvNTT converts back into the poly space in-place.
// p => p = InvNTT(p)
func (p Polynomial) InvNTT() Polynomial {
	if p.Ref.IsNTT {
		p.Ref.IsNTT = false
		p.BaseRing.InvNTT(p.Ref, p.Ref)
	}
	return p
}

// LRot performs a left rotation on the coefficients by the given positive amount in-place.
func (p Polynomial) LRot(amount int) Polynomial {
	p.InvNTT()
	p.BaseRing.Shift(p.Ref, amount, p.Ref)
	return p
}

// RRot performs a right rotation on the coefficients by the given positive amount in-place.
func (p Polynomial) RRot(amount int) Polynomial {
	leftShiftAmount := p.Ref.N() - amount
	p.LRot(leftShiftAmount)
	return p
}

// SetCoefficient updates the coefficient at the given index `i` and propagates the update through the levels.
func (p Polynomial) SetCoefficient(i int, newValue uint64) Polynomial {
	p.InvNTT()
	for lvl := 0; lvl <= p.Ref.Level(); lvl++ {
		p.Ref.Coeffs[lvl][i] = newValue
	}
	p.BaseRing.Reduce(p.Ref, p.Ref)
	return p
}

// Coeffs returns the vector representation of the zero-level coefficients of this polynomial.
func (p Polynomial) Coeffs(q *big.Int) IntVector {
	p.InvNTT()
	v := make([]RingElement, p.Ref.N())
	for i := 0; i < len(v); i++ {
		v[i] = NewModInt(int64(p.Ref.Coeffs[0][i]), q)
	}
	return NewVectorFromSlice(v).AsIntVec()
}

// Coeff returns the ith coefficient.
func (p Polynomial) Coeff(i int) uint64 {
	p.InvNTT()
	return p.Ref.Coeffs[0][i]
}

func (p Polynomial) String() string {
	p.InvNTT()
	strs := make([]string, 0, p.Ref.N())
	for _, e := range p.Ref.Coeffs[0] {
		strs = append(strs, fmt.Sprint(e))
	}
	return fmt.Sprint("Poly{" + strings.Join(strs[:5], ", ") + ", ..., " + strings.Join(strs[len(strs)-5:], ", ") + "}")
}
