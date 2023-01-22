package main

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/ldsec/codeBase/commitment/relations"
)

func getRandomKeyGenPublicParams(config relations.RelationsConfig) relations.KeyGenPublicParams {
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	}, config.Ring.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	}, config.Ring.BaseRing)
	T := fastmath.LoadNTTTransform("NTTTransform.test", config.Ring.BaseRing)
	p := config.P
	params := relations.KeyGenPublicParams{P1: p1, A1: A1, A2: A2, T: T, P: p}
	return params
}

func main() {
	config := relations.NewRelationsConfig(crypto.GetDefaultCryptoConfig())
	params := getRandomKeyGenPublicParams(config)

	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	r := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)
	k := crypto.GetAjtaiKappa(comP, comQ, params.P, config.Ring.BaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := relations.GenerateKeyGenRelation(s, r, e, k, params, config)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "rebasing")
	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebased(commitmentRing)
	// commitmentRing := config.Ring
	// rebasedRel := rel
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel).
		WithABP(128, config.Ring.Q, fastmath.NewSlice(config.Ring.D*6, config.Ring.D*7)).
		WithTernarySlice(fastmath.NewSlice(0, config.Ring.D)).
		WithReplication(4)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "public parameters creation")
	protocolParams := fastens20.GeneratePublicParameters(protocolConfig)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol execution")
	if !fastens20.Execute(rebasedRel.S, protocolParams) {
		panic("execution failed")
	}
	e0.LogExecEnd()
}
