package crypto

import (
	"math"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type RLWEConfig struct {
	fastmath.RingParams
	LogD         int
	Delta        int // Error width (max infinity norm)
	LogDelta     int
	ErrorSampler fastmath.PolySampler
}

func NewRLWEConfig(errorWidth int, ringParams fastmath.RingParams) RLWEConfig {
	logBeta := int(math.Log2(float64(errorWidth)))
	logD := int(math.Log2(float64(ringParams.D)))
	// Create the error sampler
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng for rlwe error sampler: %s")
	}
	originalGaussianSampler := ring.NewGaussianSampler(prng, ringParams.BaseRing, ringParams.Sigma, errorWidth)
	gaussianSampler := fastmath.NewAugmentedGaussianSampler(originalGaussianSampler, uint64(errorWidth), ringParams.BaseRing)
	return RLWEConfig{
		RingParams:   ringParams,
		Delta:        errorWidth,
		LogD:         logD,
		LogDelta:     logBeta,
		ErrorSampler: gaussianSampler,
	}
}

// RLWESample returns p0 and e where p0 = -p1 * s + e
func RLWESample(p1, s *fastmath.Poly, config RLWEConfig) (*fastmath.PolyNTT, *fastmath.PolyNTT) {
	e := fastmath.NewRandomPoly(config.ErrorSampler, config.BaseRing).NTT()
	return p1.Copy().Neg().NTT().Mul(s.Copy().NTT()).Add(e), e
}

// RLWEErrorDecomposition returns the ternary decomposition of the error {e_i}, b.
func RLWEErrorDecomposition(err *fastmath.Poly, config RLWEConfig) (*fastmath.IntMatrix, *fastmath.IntVec) {
	eCoeffs := err.Coeffs()
	eDecomp, ternaryBasis := fastmath.TernaryDecomposition(eCoeffs, config.LogDelta, config.BaseRing)
	return eDecomp, ternaryBasis
}

// NewIndependentRLWE construct an independent RLWE equation p0 = -p1 * s + err
func NewIndependentRLWE(p0 *fastmath.PolyNTT, p1, s, err *fastmath.Poly, T *fastmath.IntMatrix, config RLWEConfig) *LinearEquation {
	eqn := NewLinearEquation(p0.Coeffs(), T.Cols())
	Tp1 := T.Copy().(*fastmath.IntMatrix).DiagMulMat(p1.Copy().NTT().Neg().Coeffs())
	eqn.AppendTerm(Tp1, s.Coeffs()).
		AppendRLWEErrorDecompositionSum(err, T, config)
	return eqn
}

// NewIndependentRLWE construct an independent RLWE equation p0 = A * s + err
func NewIndependentRLWEWithMatrix(p0 *fastmath.PolyNTT, A *fastmath.IntMatrix, s, err *fastmath.Poly, T *fastmath.IntMatrix, config RLWEConfig) *LinearEquation {
	eqn := NewLinearEquation(p0.Coeffs(), T.Cols())
	TA := T.Copy().Hadamard(A)
	eqn.AppendTerm(TA, s.Coeffs()).
		AppendRLWEErrorDecompositionSum(err, T, config)
	return eqn
}

func (eqn *LinearEquation) AppendRLWEErrorDecompositionSum(err *fastmath.Poly, T *fastmath.IntMatrix, params RLWEConfig) *LinearEquation {
	e, b := RLWEErrorDecomposition(err, params)
	for i := 0; i < b.Size(); i++ {
		eqn.AppendTerm(T.Copy().(*fastmath.IntMatrix).ScaleCoeff(b.GetCoeff(i)), e.RowView(i))
	}
	return eqn
}
