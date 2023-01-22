package crypto

import (
	"math"
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type RLWEParameters struct {
	BaseRing *ring.Ring
	Beta     *big.Int
	LogD     int
	LogBeta  int
}

func NewRLWEParameters(q *big.Int, d int, beta uint64, baseRing *ring.Ring) RLWEParameters {
	logBeta := int(math.Log2(float64(beta)))
	logD := int(math.Log2(float64(d)))
	return RLWEParameters{
		BaseRing: baseRing,
		Beta:     big.NewInt(int64(beta)),
		LogD:     logD,
		LogBeta:  logBeta,
	}
}

// GetRLWEP0 returns -p1 * s + e in NTT domain.
func GetRLWEP0(p1, s, e *fastmath.Poly) *fastmath.PolyNTT {
	return p1.Copy().Neg().NTT().Mul(s.Copy().NTT()).Add(e.Copy().NTT())
}

// RLWEErrorDecomposition returns the ternary decomposition of the error {e_i}, b.
func RLWEErrorDecomposition(err *fastmath.Poly, params RLWEParameters) (*fastmath.IntMatrix, *fastmath.IntVec) {
	eCoeffs := err.Coeffs()
	// Decompose with the 0th level modulus.
	eDecomp, ternaryBasis := fastmath.TernaryDecomposition(eCoeffs, params.Beta, params.LogBeta, params.BaseRing)
	return eDecomp, ternaryBasis
}

// NewIndependentRLWE construct an independent RLWE equation p0 = -p1 * s + e.
func NewIndependentRLWE(p0 *fastmath.PolyNTT, p1, s, e *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearEquation {
	eqn := NewLinearEquation(p0.Coeffs(), T.Cols())
	Tp1 := T.Copy().(*fastmath.IntMatrix).DiagMulMat(p1.Copy().NTT().Neg().Coeffs())
	eqn.AppendTerm(Tp1, s.Coeffs()).
		AppendRLWEErrorDecompositionSum(e, T, params)
	return eqn
}

func (eqn *LinearEquation) AppendRLWEErrorDecompositionSum(err *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearEquation {
	e, b := RLWEErrorDecomposition(err, params)
	for i := 0; i < b.Size(); i++ {
		eqn.AppendTerm(T.Copy().(*fastmath.IntMatrix).ScaleCoeff(b.GetCoeff(i)), e.RowView(i))
	}
	return eqn
}
