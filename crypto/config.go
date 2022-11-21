package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type Config struct {
	fastmath.RingParams
	P               *big.Int // Ajtai prime mod (~20 bit prime)
	M               int      // # rows
	N               int      // # cols, must be >= D
	K               int      // Repetition rate
	Delta1          int      // Width of the uniform distribution
	Lambda          int      // M-LWE dimension
	Kappa           int      // M-SIS dimension
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
	defaultRing := fastmath.BFVZeroLevelRing()
	delta1 := 16
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	uniformSampler := ring.NewUniformSampler(prng, defaultRing.BaseRing)
	ternarySampler := ring.NewTernarySampler(prng, defaultRing.BaseRing, 1.0/3.0, false)
	gaussianSampler := ring.NewGaussianSampler(prng, defaultRing.BaseRing, defaultRing.Sigma, delta1)
	return Config{
		RingParams:      defaultRing,
		P:               big.NewInt(5857),
		N:               defaultRing.D,
		M:               16,
		K:               1,
		Delta1:          delta1,
		Lambda:          1,
		Kappa:           1,
		TernaryLength:   defaultRing.D,
		UniformSampler:  uniformSampler,
		TernarySampler:  ternarySampler,
		GaussianSampler: gaussianSampler,
	}
}
