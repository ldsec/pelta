package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type RelationsConfig struct {
	Delta1          uint64
	Beta            uint64
	P               *big.Int
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
	Ring            fastmath.RingParams
	RLWEParams      crypto.RLWEConfig
	AjtaiParams     crypto.AjtaiConfig
}

func NewRelationsConfig(cryptoConfig crypto.CryptoConfig) RelationsConfig {
	rlweParams := crypto.NewRLWEConfig(cryptoConfig.RLWEErrorWidth, cryptoConfig.RingParams)
	ajtaiParams := crypto.NewAjtaiConfig(cryptoConfig.AjtaiMod, cryptoConfig.RingParams)
	return RelationsConfig{
		Delta1:          cryptoConfig.Delta1,
		Beta:            cryptoConfig.Beta,
		P:               cryptoConfig.AjtaiMod,
		UniformSampler:  cryptoConfig.UniformSampler,
		TernarySampler:  cryptoConfig.TernarySampler,
		GaussianSampler: cryptoConfig.GaussianSampler,
		Ring:            cryptoConfig.RingParams,
		RLWEParams:      rlweParams,
		AjtaiParams:     ajtaiParams,
	}
}
