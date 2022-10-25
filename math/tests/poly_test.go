package tests

import (
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

// Creates a new polynomial with 0-level coefficients 0, 1, ..., N-1
func newTestPolynomial() (math.Polynomial, *ring.Ring) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	// Create a polynomial with coefficients 0, 1, 2, ..., N-1
	p0 := math.NewZeroPolynomial(baseRing)
	for i := 0; i < p0.Ref.N(); i++ {
		p0.SetCoefficient(i, uint64(i))
	}
	return p0, baseRing
}

// Creates a new polynomial with 0-level coefficients as the given vector.
func newTestPolynomialFrom(coeffs []uint64) (math.Polynomial, *ring.Ring) {
	// Initialize the ring.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	baseRing := ringParams.RingQP().RingQ
	p0 := math.NewZeroPolynomial(baseRing)
	for i := 0; i < len(coeffs); i++ {
		p0.SetCoefficient(i, coeffs[i])
	}
	return p0, baseRing
}

func TestNTT(t *testing.T) {
	// TODO
}

func TestInvNTT(t *testing.T) {
	// TODO
}

func TestPolyLRot(t *testing.T) {
	p0, _ := newTestPolynomial()
	maxShift := p0.Ref.N()
	for lShiftAmount := 0; lShiftAmount < maxShift; lShiftAmount++ {
		//fmt.Printf("Poly.LRot=%d\n", lShiftAmount)
		// Copy the polynomial.
		p := p0.Copy().(*math.Polynomial)
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

func TestPolyRRot(t *testing.T) {
	p0, _ := newTestPolynomial()
	maxShift := p0.Ref.N()
	for rShiftAmount := 0; rShiftAmount < maxShift; rShiftAmount++ {
		//fmt.Printf("Poly.LRot=%d\n", rShiftAmount)
		// Copy the polynomial.
		p := p0.Copy().(*math.Polynomial)
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

func TestPolyPerm(t *testing.T) {
	p0, baseRing := newTestPolynomial()
	k := 4
	sig := math.NewAutomorphism(int64(baseRing.N), int64(k))
	for exp := 1; exp < 10000; exp++ {
		p := sig.Permute(int64(exp), p0)
		// Compute (galEl ^ exp)
		newExpMult := sig.Exponent(uint64(exp))
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

func TestPolyPermInv(t *testing.T) {
	p0, baseRing := newTestPolynomial()
	k := 4
	sig := math.NewAutomorphism(int64(baseRing.N), int64(k))
	for exp := 1; exp < 10000; exp++ {
		p := sig.Permute(-int64(exp), sig.Permute(int64(exp), p0))
		// Check that p = p0
		if !p.Eq(p0) {
			t.Errorf("Poly.InvPerm: Inequality")
		}
	}
}

func TestPolyMul(t *testing.T) {
	// 4x^2 + 2x + 7
	p0, _ := newTestPolynomialFrom([]uint64{7, 2, 4})
	// 19x^4 + 2x^2 + 3x + 6
	p1, _ := newTestPolynomialFrom([]uint64{6, 3, 2, 0, 19})
	// 76x^6 + 38x^5 + 141x^4 + 16x^3 + 44x^2 + 33x + 42
	p3, _ := newTestPolynomialFrom([]uint64{42, 33, 44, 16, 141, 38, 76})
	p0.Mul(p1)
	if !p0.Eq(p3) {
		t.Errorf("Poly.Mul")
	}
}

func TestPolyPow(t *testing.T) {
	// 19x^4 + 2x^2 + 3x + 6
	p0, _ := newTestPolynomialFrom([]uint64{6, 3, 2, 0, 19})
	// 361x^8 + 76x^6 + 114x^5 + 232x^4 + 12x^3 + 33x^2 + 36x + 36
	p0_2, _ := newTestPolynomialFrom([]uint64{36, 36, 33, 12, 232, 114, 76, 0, 361})
	if !p0.Copy().Pow(2).Eq(p0_2) {
		t.Errorf("Poly.Pow")
	}
}

func TestTrace(t *testing.T) {
	// TODO
}
