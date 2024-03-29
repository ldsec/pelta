package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestPolyPerm(t *testing.T) {
	baseRing := getBaseRing()
	p0 := fastmath.NewPoly(baseRing)
	for i := 0; i < p0.N(); i++ {
		p0.SetForce(i, uint64(i+1))
	}
	k := 4
	sig := fastmath.NewAutomorphism(uint64(baseRing.N), uint64(k))
	for exp := 1; exp < 10000; exp++ {
		p := p0.Permute(int64(exp), sig)
		// Compute (galEl ^ exp)
		newExpMult := sig.Exponent(int64(exp))
		// Perform checks X^i => X^(i*galEl^exp)
		// Note that we check consistency only for the first p.Ref.N() / (galEl^exp) elements
		for i := 0; i < int(uint64(p.N())/newExpMult); i++ {
			newExp := int(uint64(i) * newExpMult % uint64(p.N()))
			if p0.GetLevel(i, 0) != (p0.GetLevel(newExp, 0) % baseRing.ModulusAtLevel[0].Uint64()) {
				t.Errorf("inconsistency at p[%d] = %d, p'[%d] = %d", i, p0.GetLevel(i, 0), newExp, p.GetLevel(newExp, 0))
			}
		}
	}
}

func TestPolyPermInv(t *testing.T) {
	baseRing := getBaseRing()
	p0 := fastmath.NewPoly(baseRing)
	for i := 0; i < p0.N(); i++ {
		p0.SetForce(i, uint64(i+1))
	}
	k := 4
	sig := fastmath.NewAutomorphism(uint64(baseRing.N), uint64(k))
	for exp := 1; exp < 10000; exp++ {
		// Permute back.
		p := p0.Permute(-int64(exp), sig)
		p = p.Permute(int64(exp), sig)
		// Check that p = p0
		if !p.Eq(p0) {
			t.Errorf("equivalence violated at exp=%d", exp)
		}
	}
}
