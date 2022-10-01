package math

import "github.com/ldsec/lattigo/v2/ring"

type Polynomial struct {
	Ref *ring.Poly
}

func NewPolynomial(baseRing *ring.Ring) Polynomial {
	return Polynomial{baseRing.NewPoly()}
}

// DeepCopy copies the polynomial.
func (p Polynomial) DeepCopy() Polynomial {
	return Polynomial{p.Ref.CopyNew()}
}

// out := op(p)
func (p Polynomial) ApplyUnaryOp(op UnaryOperator, out Polynomial, baseRing *ring.Ring) {
	op.Apply(p, out, baseRing)
}

// out := op(p, q)
func (p Polynomial) ApplyBinaryOp(op BinaryOperator, r Polynomial, out Polynomial, baseRing *ring.Ring) {
	op.Apply(p, r, out, baseRing)
}
