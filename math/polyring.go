package math

import "github.com/ldsec/lattigo/v2/ring"

type Polynomial struct {
	Ref      *ring.Poly
	BaseRing *ring.Ring
}

func NewPolynomial(baseRing *ring.Ring) Polynomial {
	return Polynomial{baseRing.NewPoly(), baseRing}
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

func (p Polynomial) MulAdd(q RingElement, out RingElement) RingElement {
	p.BaseRing.MulCoeffsAndAdd(p.Ref, q.(Polynomial).Ref, out.(Polynomial).Ref)
	return p
}

func (p Polynomial) Neg() RingElement {
	p.BaseRing.Neg(p.Ref, p.Ref)
	return p
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

// Scale scales the coefficients of the polynomial in-place.
// p => p = c*p
func (p Polynomial) Scale(factor uint64) Polynomial {
	p.BaseRing.MulScalar(p.Ref, factor, p.Ref)
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
