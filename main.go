package main

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func KeyGen(s, e, r *fastmath.Poly) {

}

func KeySwitchCollDec() {

}

func PubKeySwitch() {

}

func main() {
	// Create the main ring.
	bfvRing := fastmath.BFVZeroLevelRing()
	delta1 := 16
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	uniformSampler := ring.NewUniformSampler(prng, bfvRing.BaseRing)
	ternarySampler := ring.NewTernarySampler(prng, bfvRing.BaseRing, 1.0/3.0, false)
	// gaussianSampler := ring.NewGaussianSampler(prng, bfvRing.BaseRing, bfvRing.Sigma, delta1)
	p := big.NewInt(5857)
	beta := uint64(delta1)

	// Ajtai
	s := fastmath.NewRandomPoly(ternarySampler, bfvRing.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(bfvRing.D, bfvRing.BaseRing)
	A1 := fastmath.NewRandomIntMatrix(bfvRing.D, bfvRing.D, bfvRing.Q, bfvRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrix(bfvRing.D, bfvRing.D, bfvRing.Q, bfvRing.BaseRing)
	aj := crypto.NewAjtaiCommitment(A1, A2, s.Coeffs(), r, p, bfvRing.BaseRing)

	// RLWE
	p1 := fastmath.NewRandomPoly(uniformSampler, bfvRing.BaseRing)
	e := fastmath.NewRandomPoly(ternarySampler, bfvRing.BaseRing)
	rlweParams := crypto.NewRLWEParameters(bfvRing.Q, bfvRing.D, beta, bfvRing.BaseRing)
	rlwe := crypto.NewRLWERelation(p1, s, e, rlweParams)
	linRel := rlwe.ToLinearRelation()
	aj.EmbedIntoLinearRelation(&linRel, bfvRing.D, bfvRing.Q, bfvRing.BaseRing)
}
