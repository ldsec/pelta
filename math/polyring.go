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
	poly.SetCoefficient(0, 1)
	return poly
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
	p.BaseRing.Reduce(p.Ref, p.Ref)
	return p
}

func (p Polynomial) Sub(q RingElement) RingElement {
	return p.Add(q.Copy().Neg())
}

func (p Polynomial) Mul(q RingElement) RingElement {
	p.NTT()
	q.(Polynomial).NTT()
	p.BaseRing.MulCoeffs(p.Ref, q.(Polynomial).Ref, p.Ref)
	p.InvNTT()
	q.(Polynomial).InvNTT()
	return p
}

func (p Polynomial) MulAdd(q RingElement, out RingElement) {
	p.NTT()
	q.(Polynomial).NTT()
	out.(Polynomial).NTT()
	p.BaseRing.MulCoeffsAndAdd(p.Ref, q.(Polynomial).Ref, out.(Polynomial).Ref)
	p.InvNTT()
	q.(Polynomial).InvNTT()
	out.(Polynomial).InvNTT()
}

func (p Polynomial) Neg() RingElement {
	p.BaseRing.Neg(p.Ref, p.Ref)
	p.BaseRing.Reduce(p.Ref, p.Ref)
	return p
}

// Scale scales the coefficients of the polynomial in-place.
// p => p = c*p
func (p Polynomial) Scale(factor uint64) RingElement {
	p.NTT()
	p.BaseRing.MulScalar(p.Ref, factor, p.Ref)
	p.InvNTT()
	return p
}

func (p Polynomial) Pow(exp uint64) RingElement {
	if exp == 0 {
		return NewOnePolynomial(p.BaseRing)
	}
	out := p.Copy()
	for i := uint64(1); i < exp; i++ {
		out.Mul(p)
	}
	p.BaseRing.Reduce(out.(Polynomial).Ref, out.(Polynomial).Ref)
	p.Ref.CopyValues(out.(Polynomial).Ref)
	return p
}

func (p Polynomial) Eq(el RingElement) bool {
	for i, c1 := range p.Ref.Coeffs[0] {
		if c1 != el.(Polynomial).Ref.Coeffs[0][i] {
			return false
		}
	}
	return true
}

// Coeffs returns the vector representation of the zero-level coefficients of this polynomial.
func (p Polynomial) Coeffs(q *big.Int) IntVector {
	v := make([]RingElement, p.Ref.N())
	for i := 0; i < len(v); i++ {
		v[i] = NewModInt(int64(p.Ref.Coeffs[0][i]), q)
	}
	return NewVectorFromSlice(v).AsIntVec()
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

// Zero sets all the coefficients of this polynomial to zero.
func (p Polynomial) Zero() RingElement {
	p.Ref.Zero()
	return p
}

// One sets the first coefficient of this polynomial to one and the remaining to zero.
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

// LRot performs a left rotation on the coefficients by the given positive amount in-place.
func (p Polynomial) LRot(amount int) Polynomial {
	p.BaseRing.Shift(p.Ref, amount, p.Ref)
	return p
}

// RRot performs a right rotation on the coefficients by the given positive amount in-place.
func (p Polynomial) RRot(amount int) Polynomial {
	leftShiftAmount := p.Ref.N() - amount
	p.LRot(leftShiftAmount)
	return p
}

// Coeff returns the ith coefficient.
func (p Polynomial) Coeff(i int) uint64 {
	return p.Ref.Coeffs[0][i]
}

func (p Polynomial) String() string {
	strs := make([]string, 0, p.Ref.N())
	for _, e := range p.Ref.Coeffs[0] {
		strs = append(strs, fmt.Sprint(e))
	}
	return fmt.Sprint("Poly{" + strings.Join(strs[:5], ", ") + ", ..., " + strings.Join(strs[len(strs)-5:], ", ") + "}")
}

// Automorphism represents a Galois automorphism over polynomial rings.
type Automorphism struct {
	g int64 // The automorphism generator
	d int64 // Polynomial degree
}

func NewAutomorphism(d, k int64) Automorphism {
	return Automorphism{2*d/k + 1, d}
}

// Permute returns the permutation (sig^exp)(p) s.t. X^i => X^(i*(g^exp)) where g = (2N/k + 1)
func (aut Automorphism) Permute(exp int64, p Polynomial) Polynomial {
	var gen uint64
	if exp >= 0 {
		gen = aut.Exponent(uint64(exp))
	} else {
		// Get the additive inverse of g^exp under mod d => (exp mod d) for exp < 0
		invExp := big.NewInt(0).Mod(big.NewInt(exp), big.NewInt(aut.d)).Uint64()
		gen = aut.Exponent(invExp)
	}
	// Write the permuted result on out
	out := NewZeroPolynomial(p.BaseRing)
	p.BaseRing.Permute(p.Ref, gen, out.Ref)
	p.BaseRing.Reduce(out.Ref, out.Ref)
	return out
}

// Trace calculates the trace of this polynomial: sig^0(p) + sig^1(p) + ... + sig^k(p) and returns the result.
func (aut Automorphism) Trace(p Polynomial, k int) Polynomial {
	return NewVectorFromSize(k).Populate(
		func(v int) RingElement {
			return aut.Permute(int64(v), p)
		}).Sum().(Polynomial)
}

// Exponent computes the exponent multiplier g^exp
func (aut Automorphism) Exponent(exp uint64) uint64 {
	return big.NewInt(0).Exp(big.NewInt(aut.g), big.NewInt(int64(exp)), nil).Uint64()
}
