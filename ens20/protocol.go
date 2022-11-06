package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	"math/big"
)

// Settings represents the protocol settings.
type Settings struct {
	D               int      // deg(X^D + 1), a power of two
	LogD            int      // log(D)
	Q               *big.Int // Prime mod
	M               int      // # rows
	N               int      // # cols, must be >= D
	K               int      // Repetition rate
	Delta1          int      // Width of the uniform distribution
	Lambda          int      // M-LWE dimension
	Kappa           int      // M-SIS dimension
	BaseRing        *ring.Ring
	UniformSampler  math.PolySampler
	TernarySampler  math.PolySampler
	GaussianSampler math.PolySampler
}

// NumSplits returns n/d
func (s Settings) NumSplits() int {
	return s.N / s.D
}

// Beta returns the norm Gaussian limit.
func (s Settings) Beta() int {
	return s.Delta1
}

func GetSimpleTestSettings() Settings {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: %s")
	}
	settings := Settings{
		D:      ringParams.N(),
		LogD:   ringParams.LogN(),
		Q:      ringParams.RingQP().RingQ.ModulusAtLevel[0],
		N:      ringParams.N(),
		M:      16,
		K:      1,
		Delta1: 16,
		Lambda: 1,
		Kappa:  1,
	}
	// Initialize the ring.
	settings.BaseRing = ringParams.RingQP().RingQ
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	settings.UniformSampler = ring.NewUniformSampler(prng, settings.BaseRing)
	settings.TernarySampler = ring.NewTernarySampler(prng, settings.BaseRing, 1.0/3.0, false)
	settings.GaussianSampler = ring.NewGaussianSampler(prng, settings.BaseRing, ringParams.Sigma(), settings.Delta1)
	return settings
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	A  algebra.Matrix
	U  rings.ZIntVector
	B0 algebra.Matrix
	B  algebra.Matrix
}

func NewDummyPublicParameters(s rings.ZIntVector, settings Settings) PublicParams {
	A := math.NewRandomIntegerMatrix(settings.M, settings.N, settings.Q)
	// As = U
	u := rings.NewZIntVec(A.MulVec(s.AsVec()))
	bSize := settings.NumSplits() + 3
	B0 := math.NewRandomPolynomialMatrix(settings.Kappa, settings.Lambda+settings.Kappa+bSize, settings.BaseRing, settings.UniformSampler)
	b := math.NewRandomPolynomialMatrix(bSize, B0.Cols(), settings.BaseRing, settings.UniformSampler)

	return PublicParams{
		A:  A,
		U:  u,
		B0: B0,
		B:  b,
	}
}
