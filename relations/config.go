package relations

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

var mainRing = fastmath.BFVFullRingCustom(bfv.PN13QP218, 13)
var rebaseRing = fastmath.BFVFullRingCustom(bfv.PN13QP218, 7)

func GetDefaultAjtaiConfig() crypto.AjtaiConfig {
	p := big.NewInt(int64(uint64(1<<20 - 3)))
	return crypto.NewAjtaiConfig(p, mainRing)
}

func GetDefaultRLWEConfig() crypto.RLWEConfig {
	rlweErrorWidth := 16
	return crypto.NewRLWEConfig(rlweErrorWidth, mainRing)
}

func GetDefaultProtocolConfig(m, n int) fastens20.ProtocolConfig {
	delta1 := uint64(1 << 24)
	return fastens20.NewProtocolConfig(rebaseRing, delta1, m, n).
		WithReplication(4).
		WithSecurityParameters(9, 10)
}
