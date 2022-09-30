package math

import "github.com/ldsec/lattigo/v2/ring"

// BinaryOperator represents a binary operation between two polynomials.
type BinaryOperator interface {
	Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring)
}

// UnaryOperator represents a unary operation on a polynomial.
type UnaryOperator interface {
	Apply(a *ring.Poly, out *ring.Poly, baseRing *ring.Ring)
}

// NTTMul represents a multiplication of two polynomials that are in the NTT domain. The output is also in NTT-domain.
// Represents out = a * b
type NTTMul struct{}

func (_ NTTMul) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.MulCoeffs(a, b, out)
}

// PolyMul represents a multiplication of two polynomials that are in poly domain. The output is also in poly domain.
// Represents out = a * b
type PolyMul struct{}

func (_ PolyMul) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.NTT(a, a)
	baseRing.NTT(b, b)
	NTTMul{}.Apply(a, b, out, baseRing)
	baseRing.InvNTT(a, a)
	baseRing.InvNTT(b, b)
	baseRing.InvNTT(out, out)
}

// NTTMulAdd represents a multiplication of two polynomials that are in the NTT domain. The output is accumulated in
// the `out`. The output is also in NTT-domain.
// Represents out += a * b.
type NTTMulAdd struct{}

func (_ NTTMulAdd) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.MulCoeffsAndAdd(a, b, out)
}

// PolyMulAdd represents a multiplication of two polynomials that are in the NTT domain. The output is accumulated in
// the `out`. The output is also in NTT-domain.
// Represents out += a * b.
type PolyMulAdd struct{}

func (_ PolyMulAdd) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.NTT(a, a)
	baseRing.NTT(b, b)
	NTTMulAdd{}.Apply(a, b, out, baseRing)
	baseRing.InvNTT(a, a)
	baseRing.InvNTT(b, b)
	baseRing.InvNTT(out, out)
}

// SimpleAdd represents a coefficient wise addition of two polynomials.
type SimpleAdd struct{}

func (_ SimpleAdd) Apply(a *ring.Poly, b *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.Add(a, b, out)
}

// NTT converts a polynomial to its NTT form.
type NTT struct{}

func (_ NTT) Apply(a *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.NTT(a, out)
}

// InvNTT converts a polynomial from its NTT form to its normal form.
type InvNTT struct{}

func (_ InvNTT) Apply(a *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.InvNTT(a, out)
}

// Scale multiplies the coefficients of a polynomial with a given factor.
type Scale struct {
	Factor uint64
}

func (s Scale) Apply(a *ring.Poly, out *ring.Poly, baseRing *ring.Ring) {
	baseRing.MulScalar(a, s.Factor, out)
}
