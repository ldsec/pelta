package main

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

var mainRing fastmath.RingParams = fastmath.BFVFullRing()
var rebaseRing fastmath.RingParams = fastmath.BFVFullShortCommtRing(7)

func GetDefaultAjtaiConfig() crypto.AjtaiConfig {
	p := big.NewInt(5857)
	return crypto.NewAjtaiConfig(p, mainRing)
}

func GetDefaultRLWEConfig() crypto.RLWEConfig {
	rlweErrorWidth := 128
	return crypto.NewRLWEConfig(rlweErrorWidth, mainRing)
}

func GetDefaultProtocolConfig(m, n int) fastens20.ProtocolConfig {
	delta1 := uint64(1 << 24)
	return fastens20.NewProtocolConfig(rebaseRing, delta1, m, n)
}
