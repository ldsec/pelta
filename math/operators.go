package math

import "github.com/ldsec/lattigo/v2/ring"

// BinaryOperator represents a binary operation between two polynomials.
type BinaryOperator interface {
	Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring)
}

// UnaryOperator represents a unary operation on a polynomial.
type UnaryOperator interface {
	Apply(a Polynomial, out Polynomial, baseRing *ring.Ring)
}

// NTTMul represents a multiplication of two polynomials that are in the NTT domain. The output is also in NTT-domain.
// Represents out = a * b
type NTTMul struct{}

func (_ NTTMul) Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.MulCoeffs(a.Ref, b.Ref, out.Ref)
}

// PolyMul represents a multiplication of two polynomials that are in poly domain. The output is also in poly domain.
// Represents out = a * b
type PolyMul struct{}

func (_ PolyMul) Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring) {
	NTT{}.Apply(a, a, baseRing)
	NTT{}.Apply(b, b, baseRing)
	NTTMul{}.Apply(a, b, out, baseRing)
	InvNTT{}.Apply(a, a, baseRing)
	InvNTT{}.Apply(b, b, baseRing)
	InvNTT{}.Apply(out, out, baseRing)
}

// NTTMulAdd represents a multiplication of two polynomials that are in the NTT domain. The output is accumulated in
// the `out`. The output is also in NTT-domain.
// Represents out += a * b.
type NTTMulAdd struct{}

func (_ NTTMulAdd) Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.MulCoeffsAndAdd(a.Ref, b.Ref, out.Ref)
}

// PolyMulAdd represents a multiplication of two polynomials that are in the NTT domain. The output is accumulated in
// the `out`. The output is also in NTT-domain.
// Represents out += a * b.
type PolyMulAdd struct{}

func (_ PolyMulAdd) Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring) {
	NTT{}.Apply(a, a, baseRing)
	NTT{}.Apply(b, b, baseRing)
	NTTMulAdd{}.Apply(a, b, out, baseRing)
	InvNTT{}.Apply(a, a, baseRing)
	InvNTT{}.Apply(b, b, baseRing)
	InvNTT{}.Apply(out, out, baseRing)
}

// SimpleAdd represents a coefficient wise addition of two polynomials.
type SimpleAdd struct{}

func (_ SimpleAdd) Apply(a Polynomial, b Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.Add(a.Ref, b.Ref, out.Ref)
}

// NTT converts a polynomial to its NTT form.
type NTT struct{}

func (_ NTT) Apply(a Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.NTT(a.Ref, out.Ref)
}

// InvNTT converts a polynomial from its NTT form to its normal form.
type InvNTT struct{}

func (_ InvNTT) Apply(a Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.InvNTT(a.Ref, out.Ref)
}

// Scale multiplies the coefficients of a polynomial with a given factor.
type Scale struct {
	Factor uint64
}

func (s Scale) Apply(a Polynomial, out Polynomial, baseRing *ring.Ring) {
	baseRing.MulScalar(a.Ref, s.Factor, out.Ref)
}
