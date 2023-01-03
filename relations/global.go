package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type GlobalConfig struct {
	Delta1          int
	Beta            uint64
	P               *big.Int
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
	RLWEParams      crypto.RLWEParameters
	Ring            fastmath.RingParams
}

func NewGlobalConfig(cryptoConfig crypto.CryptoConfig) GlobalConfig {
	beta := uint64(16)
	rlweParams := crypto.NewRLWEParameters(cryptoConfig.Q, cryptoConfig.D, beta, cryptoConfig.BaseRing)
	return GlobalConfig{
		Delta1:          cryptoConfig.Delta1,
		Beta:            beta,
		P:               cryptoConfig.P,
		UniformSampler:  cryptoConfig.UniformSampler,
		TernarySampler:  cryptoConfig.TernarySampler,
		GaussianSampler: cryptoConfig.GaussianSampler,
		RLWEParams:      rlweParams,
		Ring:            cryptoConfig.RingParams,
	}
}
