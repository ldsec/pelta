package relations

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

var MainRing = fastmath.BFVFullRingCustom(bfv.PN13QP218, 13)
var RebaseRing = fastmath.BFVFullRingCustom(bfv.PN13QP218, 7)
var Delta1 = uint64(1 << 25)
var Kappa = 9
var Lambda = 10

func GetDefaultAjtaiConfig() crypto.AjtaiConfig {
	p := big.NewInt(int64(uint64(1<<20 - 3)))
	return crypto.NewAjtaiConfig(p, MainRing)
}

func GetDefaultRLWEConfig() crypto.RLWEConfig {
	rlweErrorWidth := 16
	return crypto.NewRLWEConfig(rlweErrorWidth, MainRing)
}

func GetDefaultProtocolConfig(m, n int) fastens20.ProtocolConfig {
	conf := fastens20.NewProtocolConfig(RebaseRing, Delta1, m, n).
		WithReplication(4).
		WithSecurityParameters(Kappa, Lambda)
	//conf.T = uint64(1 << 20)
	return conf
}
