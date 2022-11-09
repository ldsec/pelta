package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Automorphism represents a Galois automorphism over polynomial rings.
type Automorphism struct {
	G uint64 // The automorphism generator
	D uint64 // Polynomial degree
}

func NewAutomorphism(d, k uint64) Automorphism {
	return Automorphism{2*d/k + 1, d}
}

// Permute returns the permutation of this polynomial, i.e, sig^exp(p).
func (p *Poly) Permute(exp int64, sig Automorphism) Poly {
	var gen uint64
	if exp >= 0 {
		gen = sig.Exponent(uint64(exp))
	} else {
		// Get the additive inverse of g^exp under mod d => (exp mod d) for exp < 0
		// TODO: optimize
		invExp := big.NewInt(0).Mod(big.NewInt(exp), big.NewInt(int64(sig.D))).Uint64()
		gen = sig.Exponent(invExp)
	}
	// Write the permuted result on out
	out := NewZeroPoly(p.baseRing)
	p.baseRing.Permute(p.ref, gen, out.ref)
	p.baseRing.Reduce(out.ref, out.ref)
	return out
}

// Trace calculates the trace of this function: sig^0(f(0)) + sig^1(f(1)) + ... + sig^k(f(k)) and returns the result.
func Trace(sig Automorphism, f func(int) Poly, k int, baseRing *ring.Ring) Poly {
	out := NewPolyVec(k, baseRing)
	out.Populate(func(v int) Poly {
		val := f(v)
		return val.Permute(int64(v), sig)
	})
	return out.Sum()
}

// Exponent computes the exponent multiplier g^exp
func (sig Automorphism) Exponent(exp uint64) uint64 {
	// TODO: optimize
	return big.NewInt(0).Exp(big.NewInt(int64(sig.G)), big.NewInt(int64(exp)), nil).Uint64()
}
