package main

import "github.com/ldsec/lattigo/v2/ring"

// BinaryOperator represents a binary operation between two polynomials.
type BinaryOperator interface {
	Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring)
}

// NTTMult represents a multiplication of two polynomials that are in the NTT domain.
// The output is also in NTT-domain.
type NTTMult struct{}

func (N NTTMult) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.MulCoeffs(a, b, out)
}

// PolyMult represents a multiplication of two polynomials that are in poly domain.
// The output is also in poly domain.
type PolyMult struct{}

func (p PolyMult) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.NTT(a, a)
	baseRing.NTT(b, b)
	NTTMult{}.Apply(a, b, out, baseRing)
	baseRing.InvNTT(a, a)
	baseRing.InvNTT(b, b)
	baseRing.InvNTT(out, out)
}
