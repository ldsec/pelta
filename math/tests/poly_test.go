package tests

import (
	"fmt"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"testing"
)

// Performs equality check between two coefficient vectors.
func coeffsAtLevelEqual(coeff1 []uint64, coeff2 []uint64) bool {
	for i := 0; i < len(coeff1); i++ {
		if coeff1[i] != coeff2[i] {
			return false
		}
	}
	return true
}

// Prints the coefficients at the level.
func printCoeffsAtLevel(coeff1 []uint64) {
	for i := 0; i < len(coeff1); i++ {
		fmt.Printf("[%d]", coeff1[i])
	}
	fmt.Println()
}

// Creates a new polynomial with coefficients 0, 1, ..., N-1
func newTestPolynomial(baseRing *ring.Ring) math.Polynomial {
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := math.NewZeroPolynomial(baseRing)
	for i := 0; i < p0.Ref.N(); i++ {
		p0.SetCoefficient(i, uint64(i))
	}
	return p0
}

func TestNTT(t *testing.T) {
	// TODO
}

func TestInvNTT(t *testing.T) {
	// TODO
}

func TestLRot(t *testing.T) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := newTestPolynomial(baseRing)
	maxShift := p0.Ref.N()
	for lShiftAmount := 0; lShiftAmount < maxShift; lShiftAmount++ {
		//fmt.Printf("Poly.LRot=%d\n", lShiftAmount)
		// Copy the polynomial.
		p := p0.Copy().(math.Polynomial)
		p.LRot(lShiftAmount)
		// Create the expected coefficients
		shiftedCoeffs := make([]uint64, p.Ref.N())
		for i := 0; i < len(shiftedCoeffs); i++ {
			shiftedCoeffs[i] = uint64((i + lShiftAmount) % p.Ref.N())
		}
		//printCoeffsAtLevel(p.Ref.Coeffs[0])
		//printCoeffsAtLevel(shiftedCoeffs)
		if !coeffsAtLevelEqual(p.Ref.Coeffs[0], shiftedCoeffs) {
			t.Errorf("Poly.LRot")
		}
	}

}

func TestRRot(t *testing.T) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := newTestPolynomial(baseRing)
	maxShift := p0.Ref.N()
	for rShiftAmount := 0; rShiftAmount < maxShift; rShiftAmount++ {
		//fmt.Printf("Poly.LRot=%d\n", rShiftAmount)
		// Copy the polynomial.
		p := p0.Copy().(math.Polynomial)
		p.RRot(rShiftAmount)
		// Create the expected coefficients
		shiftedCoeffs := make([]uint64, p.Ref.N())
		for i := 0; i < len(shiftedCoeffs); i++ {
			shiftedCoeffs[i] = uint64((p.Ref.N() - rShiftAmount + i) % p.Ref.N())
		}
		//printCoeffsAtLevel(p.Ref.Coeffs[0])
		//printCoeffsAtLevel(shiftedCoeffs)
		if !coeffsAtLevelEqual(p.Ref.Coeffs[0], shiftedCoeffs) {
			t.Errorf("Poly.LRot")
		}
	}
}

func TestPerm(t *testing.T) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	//q := baseRing.ModulusAtLevel[0]
	k := 4
	sig := math.NewAutomorphism(int64(baseRing.N), int64(k))
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := newTestPolynomial(baseRing)
	for exp := 1; exp < 10000; exp++ {
		p := sig.Permute(int64(exp), p0)
		// Compute (galEl ^ exp)
		newExpMult := sig.Exponent(uint64(exp))
		//fmt.Println("Poly.Perm: galEl^exp =", newExpMult)
		// Perform checks X^i => X^(i*galEl^exp)
		// Note that we check consistency only for the first p.Ref.N() / (galEl^exp) elements
		for i := uint64(0); i < uint64(p.Ref.N())/newExpMult; i++ {
			newExp := (i * newExpMult) % uint64(p.Ref.N())
			if p0.Ref.Coeffs[0][i] != (p.Ref.Coeffs[0][newExp] % baseRing.ModulusAtLevel[0].Uint64()) {
				t.Errorf("Poly.Perm: Inconsistency at p[%d] = %d, p'[%d] = %d", i, p0.Ref.Coeffs[0][i], newExp, p.Ref.Coeffs[0][newExp])
			}
		}
	}
}

func TestPermInv(t *testing.T) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	//q := baseRing.ModulusAtLevel[0]
	k := 4
	sig := math.NewAutomorphism(int64(baseRing.N), int64(k))
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := newTestPolynomial(baseRing)
	for exp := 1; exp < 10000; exp++ {
		p := sig.Permute(-int64(exp), sig.Permute(int64(exp), p0))
		// Check that p = p0
		for i := 0; i < p.Ref.N(); i++ {
			if p.Coeff(i) != p0.Coeff(i) {
				t.Errorf("Poly.Perm: Inconsistency at p[%d] = %d, p'[%d] = %d", i, p0.Ref.Coeffs[0][i], i, p.Ref.Coeffs[0][i])
			}
		}
	}
}

func TestTrace(t *testing.T) {
	// TODO
}
