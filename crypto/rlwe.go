package crypto

import (
	"math"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type RLWEParameters struct {
	fastmath.RingParams
	LogD     int
	Delta    int // Error width (max infinity norm)
	LogDelta int
}

func NewRLWEParameters(delta int, ringParams fastmath.RingParams) RLWEParameters {
	logBeta := int(math.Log2(float64(delta)))
	logD := int(math.Log2(float64(ringParams.D)))
	return RLWEParameters{
		RingParams: ringParams,
		Delta:      delta,
		LogD:       logD,
		LogDelta:   logBeta,
	}
}

func GetRLWEErrorSampler(params RLWEParameters) fastmath.PolySampler {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng for rlwe error sampler: %s")
	}
	originalGaussianSampler := ring.NewGaussianSampler(prng, params.BaseRing, params.Sigma, params.Delta)
	gaussianSampler := fastmath.NewAugmentedGaussianSampler(originalGaussianSampler, uint64(params.Delta), params.BaseRing)
	return gaussianSampler
}

// RLWESample returns -p1 * s + err in NTT domain.
func RLWESample(p1, s, err *fastmath.Poly) *fastmath.PolyNTT {
	return p1.Copy().Neg().NTT().Mul(s.Copy().NTT()).Add(err.Copy().NTT())
}

// RLWEErrorDecomposition returns the ternary decomposition of the error {e_i}, b.
func RLWEErrorDecomposition(err *fastmath.Poly, params RLWEParameters) (*fastmath.IntMatrix, *fastmath.IntVec) {
	eCoeffs := err.Coeffs()
	eDecomp, ternaryBasis := fastmath.TernaryDecomposition(eCoeffs, params.LogDelta, params.BaseRing)
	return eDecomp, ternaryBasis
}

// NewIndependentRLWE construct an independent RLWE equation p0 = -p1 * s + err
func NewIndependentRLWE(p0 *fastmath.PolyNTT, p1, s, err *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearEquation {
	eqn := NewLinearEquation(p0.Coeffs(), T.Cols())
	Tp1 := T.Copy().(*fastmath.IntMatrix).DiagMulMat(p1.Copy().NTT().Neg().Coeffs())
	eqn.AppendTerm(Tp1, s.Coeffs()).
		AppendRLWEErrorDecompositionSum(err, T, params)
	return eqn
}

func (eqn *LinearEquation) AppendRLWEErrorDecompositionSum(err *fastmath.Poly, T *fastmath.IntMatrix, params RLWEParameters) *LinearEquation {
	e, b := RLWEErrorDecomposition(err, params)
	for i := 0; i < b.Size(); i++ {
		eqn.AppendTerm(T.Copy().(*fastmath.IntMatrix).ScaleCoeff(b.GetCoeff(i)), e.RowView(i))
	}
	return eqn
}
