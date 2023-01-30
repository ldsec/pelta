package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

type CryptoConfig struct {
	fastmath.RingParams
	Delta1          uint64 // Width of the uniform distribution
	Beta            uint64 // Mask bound
	Lambda          int    // M-LWE dimension
	Kappa           int    // M-SIS dimension
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler

	AjtaiMod       *big.Int // Ajtai prime mod (~20 bit prime)
	RLWEErrorWidth int
}

func GetDefaultCryptoConfig() CryptoConfig {
	// Initialize the ring parameters.
	defaultRing := fastmath.BFVFullRing()
	delta1 := uint64(128)
	beta := delta1
	uniformSampler, ternarySampler, gaussianSampler := fastmath.GetSamplers(defaultRing, delta1)
	return CryptoConfig{
		RingParams:      defaultRing,
		AjtaiMod:        big.NewInt(5857),
		Delta1:          delta1,
		Beta:            beta,
		Lambda:          1,
		Kappa:           1,
		UniformSampler:  uniformSampler,
		TernarySampler:  ternarySampler,
		GaussianSampler: gaussianSampler,
		RLWEErrorWidth:  128,
	}
}

func GetPeltaCryptoConfig() CryptoConfig {
	// Initialize the ring parameters.
	defaultRing := fastmath.BFVFullRing()
	delta1 := uint64(1 << 24)
	beta := delta1
	p := big.NewInt(5857)
	uniformSampler, ternarySampler, gaussianSampler := fastmath.GetSamplers(defaultRing, delta1)
	return CryptoConfig{
		RingParams:      defaultRing,
		AjtaiMod:        p,
		Delta1:          delta1,
		Beta:            beta,
		Lambda:          10,
		Kappa:           9,
		UniformSampler:  uniformSampler,
		TernarySampler:  ternarySampler,
		GaussianSampler: gaussianSampler,
		RLWEErrorWidth:  128,
	}
}