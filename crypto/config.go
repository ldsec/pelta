package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type CryptoConfig struct {
	fastmath.RingParams
	P               *big.Int // Ajtai prime mod (~20 bit prime)
	Delta1          int      // Width of the uniform distribution
	Beta            int
	Lambda          int // M-LWE dimension
	Kappa           int // M-SIS dimension
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
}

func GetDefaultCryptoConfig() CryptoConfig {
	// Initialize the ring parameters.
	defaultRing := fastmath.BFVFullRing()
	delta1 := 128
	beta := delta1
	uniformSampler, ternarySampler, gaussianSampler := GetSamplers(defaultRing, delta1)
	return CryptoConfig{
		RingParams:      defaultRing,
		P:               big.NewInt(5857),
		Delta1:          delta1,
		Beta:            beta,
		Lambda:          1,
		Kappa:           1,
		UniformSampler:  uniformSampler,
		TernarySampler:  ternarySampler,
		GaussianSampler: gaussianSampler,
	}
}

// GetSamplers constructs and returns the uniform sampler, ternary sampler, and gaussian sampler
func GetSamplers(samplerRing fastmath.RingParams, delta1 int) (fastmath.PolySampler, fastmath.PolySampler, fastmath.PolySampler) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	uniformSampler := ring.NewUniformSampler(prng, samplerRing.BaseRing)
	originalTernarySampler := ring.NewTernarySampler(prng, samplerRing.BaseRing, 1.0/3.0, false)
	ternarySampler := fastmath.NewAugmentedTernarySampler(originalTernarySampler, samplerRing.BaseRing)
	gaussianSampler := ring.NewGaussianSampler(prng, samplerRing.BaseRing, samplerRing.Sigma, delta1)

	return uniformSampler, ternarySampler, gaussianSampler
}
