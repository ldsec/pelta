package crypto

import (
	"math"
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type RLWEParameters struct {
	BaseRing *ring.Ring
	Q        *big.Int
	Beta     uint64
	LogD     int
	LogBeta  int
}

func NewRLWEParameters(q *big.Int, d int, beta uint64, baseRing *ring.Ring) RLWEParameters {
	logBeta := int(math.Log2(float64(beta)))
	logD := int(math.Log2(float64(d)))
	return RLWEParameters{
		BaseRing: baseRing,
		Q:        q,
		Beta:     beta,
		LogD:     logD,
		LogBeta:  logBeta,
	}
}

// AppendIndependentRLWE appends an independent RLWE equation p0 = -p1 * s + e.
func (lrb *LinearRelationBuilder) AppendIndependentRLWE(p0 *fastmath.PolyNTT, p1, s, e *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearRelationBuilder {
	eqn := NewLinearEquation(p0.Coeffs(), T.Cols()).
		AppendTerm(T.Copy().DiagMulMat(p1.Copy().NTT().Neg().Coeffs()), s.Coeffs()).
		AppendErrorDecompositionSum(e, T, params)
	lrb.AppendEqn(eqn)
	return lrb
}

func (eqn *LinearEquation) AppendErrorDecompositionSum(err *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearEquation {
	e, b := ErrorDecomposition(err, params)
	for i := 0; i < b.Size(); i++ {
		eqn.AppendTerm(T.Copy().Scale(b.Get(i)), e.RowView(i))
	}
	return eqn
}

// ErrorDecomposition returns the ternary decomposition of the error {e_i}, b.
func ErrorDecomposition(err *fastmath.Poly, params RLWEParameters) (*fastmath.IntMatrix, *fastmath.IntVec) {
	eCoeffs := err.Coeffs()
	eDecomp, ternaryBasis := fastmath.TernaryDecomposition(eCoeffs, params.Beta, params.LogBeta, params.Q, params.BaseRing)
	return eDecomp.Transposed(), ternaryBasis
}
