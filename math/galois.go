package math

import "math/big"

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

// Trace calculates the trace of this function: sig^0(f(0)) + sig^1(f(1)) + ... + sig^k(f(k)) and returns the result.
func (aut Automorphism) Trace(f func(int) Polynomial, k int) Polynomial {
	return NewVectorFromSize(k).Populate(
		func(v int) RingElement {
			return aut.Permute(int64(v), f(v))
		}).Sum().(Polynomial)
}

// Exponent computes the exponent multiplier g^exp
func (aut Automorphism) Exponent(exp uint64) uint64 {
	return big.NewInt(0).Exp(big.NewInt(aut.g), big.NewInt(int64(exp)), nil).Uint64()
}
