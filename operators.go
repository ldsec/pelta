package main

import "github.com/ldsec/lattigo/v2/ring"

// UnaryOperator represents an unary operation on a polynomial.
type UnaryOperator interface {
	Apply(a *ring.Poly, baseRing *ring.Ring) *ring.Poly
}

// ToNTT converts a polynomial from poly domain to NTT domain.
type ToNTT struct{}

func (t ToNTT) Apply(a *ring.Poly, baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	baseRing.NTT(a, out)
	return out
}

// ToPoly converts a polynomial from NTT domain to poly domain.
type ToPoly struct{}

func (t ToPoly) Apply(a *ring.Poly, baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	baseRing.InvNTT(a, out)
	return out
}

// BinaryOperator represents a binary operation between two polynomials.
type BinaryOperator interface {
	Apply(a *ring.Poly, b *ring.Poly, baseRing *ring.Ring) *ring.Poly
}

// NTTMult represents a multiplication of two polynomials that are in the NTT domain.
// The output is also in NTT-domain.
type NTTMult struct{}

func (N NTTMult) Apply(a *ring.Poly, b *ring.Poly, baseRing *ring.Ring) *ring.Poly {
	tmp := baseRing.NewPoly()
	baseRing.MulCoeffs(a, b, tmp)
	return tmp
}

// PolyMult represents a multiplication of two polynomials that are in poly domain.
// The output is also in poly domain.
type PolyMult struct{}

func (p PolyMult) Apply(a *ring.Poly, b *ring.Poly, baseRing *ring.Ring) *ring.Poly {
	baseRing.NTT(a, a)
	baseRing.NTT(b, b)
	tmp := NTTMult{}.Apply(a, b, baseRing)
	baseRing.InvNTT(a, a)
	baseRing.InvNTT(b, b)
	baseRing.InvNTT(tmp, tmp)
	return tmp
}
