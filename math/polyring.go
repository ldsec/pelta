package math

import (
	"github.com/tuneinsight/lattigo/v4/ring"
)

type Polynomial struct {
	Ref      *ring.Poly
	BaseRing *ring.Ring
}

func NewZeroPolynomial(baseRing *ring.Ring) Polynomial {
	return Polynomial{baseRing.NewPoly(), baseRing}
}

func NewPolynomial(ref *ring.Poly, baseRing *ring.Ring) Polynomial {
	return Polynomial{ref, baseRing}
}

// Copy copies the polynomial.
func (p Polynomial) Copy() RingElement {
	return Polynomial{p.Ref.CopyNew(), p.BaseRing}
}

func (p Polynomial) Add(q RingElement) RingElement {
	p.BaseRing.Add(p.Ref, q.(Polynomial).Ref, p.Ref)
	return p
}

func (p Polynomial) Mul(q RingElement) RingElement {
	p.BaseRing.MulCoeffs(p.Ref, q.(Polynomial).Ref, p.Ref)
	return p
}

func (p Polynomial) MulAdd(q RingElement, out RingElement) {
	p.BaseRing.MulCoeffsAndAdd(p.Ref, q.(Polynomial).Ref, out.(Polynomial).Ref)
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

func (p Polynomial) Pow(exp int64) RingElement {
	out := p.Copy().One()
	expAbs := exp
	if expAbs < 0 {
		expAbs *= -1
	}
	for i := int64(0); i < expAbs; i++ {
		out.Mul(p)
	}
	// TODO Perform inversion if exp < 0
	p.Ref.CopyValues(out.(Polynomial).Ref)
	return p
}

func (p Polynomial) Eq(el RingElement) bool {
	return p.BaseRing.Equal(p.Ref, el.(Polynomial).Ref)
}

// NTT converts into the NTT space in-place.
// p => p = NTT(p)
func (p Polynomial) NTT() Polynomial {
	p.BaseRing.NTT(p.Ref, p.Ref)
	return p
}

// InvNTT converts back into the poly space in-place.
// p => p = InvNTT(p)
func (p Polynomial) InvNTT() Polynomial {
	p.BaseRing.InvNTT(p.Ref, p.Ref)
	return p
}

func (p Polynomial) Zero() RingElement {
	p.Ref.Zero()
	return p
}

func (p Polynomial) One() RingElement {
	p.Ref.Zero()
	p.Ref.Coeffs[0][0] = 1
	return p
}

// SetCoefficient updates the coefficient at the given index `i` and propagates the update through the levels.
func (p Polynomial) SetCoefficient(i int, newValue uint64) {
	for lvl := 0; lvl <= p.Ref.Level(); lvl++ {
		p.Ref.Coeffs[lvl][i] = newValue
	}
	p.BaseRing.Reduce(p.Ref, p.Ref)
}

// Perm performs a permutation sig^exp(p) s.t. X^i => X^(i*(galEl^exp)) in-place.
func (p Polynomial) Perm(galEl *ModInt, exp int64) Polynomial {
	p.BaseRing.Permute(p.Ref, galEl.Copy().Pow(exp).(*ModInt).Uint64(), p.Ref)
	return p
}

// Trace calculates the trace of this polynomial: sig^0(p) + sig^1(p) + ... + sig^k(p) and returns the result.
func (p Polynomial) Trace(galEl *ModInt, k int) Polynomial {
	return NewVectorFromSize(k).Populate(
		func(v int) RingElement {
			return p.Copy().(Polynomial).Perm(galEl, int64(v))
		}).Sum().(Polynomial)
}

// LShift performs a left shift on the coefficients by the given amount in-place.
func (p Polynomial) LShift(amount int) Polynomial {
	p.BaseRing.Shift(p.Ref, amount, p.Ref)
	return p
}

// RShift performs a right shift on the coefficients by the given amount in-place.
func (p Polynomial) RShift(amount int) Polynomial {
	leftShiftAmount := p.Ref.N() - amount
	p.LShift(leftShiftAmount)
	return p
}

// Coeff returns the ith coefficient.
func (p Polynomial) Coeff(i int) uint64 {
	return p.Ref.Coeffs[0][i]
}
