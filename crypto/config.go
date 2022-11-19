package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type Config struct {
	D               int      // deg(X^D + 1), a power of two
	LogD            int      // log(D)
	Q               *big.Int // Prime mod
	P               *big.Int // Ajtai prime mod (~20 bit prime)
	M               int      // # rows
	N               int      // # cols, must be >= D
	K               int      // Repetition rate
	Delta1          int      // Width of the uniform distribution
	Lambda          int      // M-LWE dimension
	Kappa           int      // M-SIS dimension
	BaseRing        *ring.Ring
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler

	TernaryLength int // Length of the ternary prefix of the M-SIS secret
}

// NumSplits returns n/d
func (c Config) NumSplits() int {
	return c.N / c.D
}

// Beta returns the norm Gaussian limit.
func (c Config) Beta() int {
	return c.Delta1
}

// NumTernarySplits returns the number of splits that are in ternary.
func (c Config) NumTernarySplits() int {
	return c.TernaryLength / c.D
}

func GetDefaultConfig() Config {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: %s")
	}
	settings := Config{
		D:             ringParams.N(),
		LogD:          ringParams.LogN(),
		Q:             ringParams.RingQP().RingQ.ModulusAtLevel[0],
		P:             big.NewInt(5857),
		N:             ringParams.N(),
		M:             16,
		K:             1,
		Delta1:        16,
		Lambda:        1,
		Kappa:         1,
		TernaryLength: ringParams.N(),
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
