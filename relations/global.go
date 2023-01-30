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
	RLWEParams      crypto.RLWEParameters
	AjtaiParams     crypto.AjtaiParameters
}

func NewRelationsConfig(cryptoConfig crypto.CryptoConfig) RelationsConfig {
	rlweParams := crypto.NewRLWEParameters(cryptoConfig.RLWEErrorWidth, cryptoConfig.RingParams)
	ajtaiParams := crypto.NewAjtaiParameters(cryptoConfig.AjtaiMod, cryptoConfig.RingParams)
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
