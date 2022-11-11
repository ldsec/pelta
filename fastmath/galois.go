package fastmath

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
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
func (p *Poly) Permute(exp int64, sig Automorphism) *Poly {
	var gen uint64
	if exp >= 0 {
		gen = sig.Exponent(uint64(exp))
	} else {
		// Get the additive inverse of g^exp under mod d => (exp mod d) for exp < 0
		//invExp := ring.ModExp(uint64(int64(sig.D)+exp), 1, sig.D)
		invExp := big.NewInt(0).Mod(big.NewInt(exp), big.NewInt(int64(sig.D))).Uint64()
		gen = sig.Exponent(invExp)
	}
	// Write the permuted result on out
	out := NewZeroPoly(p.baseRing)
	p.baseRing.Permute(p.ref, gen, out.ref)
	return out
}

// Trace calculates the trace of this function: sig^0(f(0)) + sig^1(f(1)) + ... + sig^{k-1}(f(k-1)) and returns the result.
func Trace(sig Automorphism, f func(int) Poly, k int, baseRing *ring.Ring) *Poly {
	out := NewPolyVec(k, baseRing)
	out.Populate(func(v int) Poly {
		fOut := f(v)
		return *fOut.Permute(int64(v), sig)
	})
	return out.Sum()
}

// Exponent computes the exponent multiplier g^exp
func (sig Automorphism) Exponent(exp uint64) uint64 {
	return big.NewInt(0).Exp(big.NewInt(int64(sig.G)), big.NewInt(int64(exp)), nil).Uint64()
	// NOTE: 2d works up to some degree. May need to tune a bit.
	//return ring.ModExp(sig.G, exp, 2*sig.D)
}
