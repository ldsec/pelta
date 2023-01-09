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
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	uniformSampler := ring.NewUniformSampler(prng, defaultRing.BaseRing)
	originalTernarySampler := ring.NewTernarySampler(prng, defaultRing.BaseRing, 1.0/3.0, false)
	ternarySampler := fastmath.NewAugmentedTernarySampler(originalTernarySampler, defaultRing.BaseRing)
	gaussianSampler := ring.NewGaussianSampler(prng, defaultRing.BaseRing, defaultRing.Sigma, delta1)
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
