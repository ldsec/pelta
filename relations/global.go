package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

type GlobalConfig struct {
	BfvRing         fastmath.RingParams
	Delta1          int
	Beta            uint64
	P               *big.Int
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
}
