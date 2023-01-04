package test

import (
	"math"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/relations"
)

func TestKeyGen(t *testing.T) {
	t.Logf("reading config")
	config := relations.NewGlobalConfig(crypto.GetDefaultCryptoConfig())
	q := config.Ring.Q
	logD := int(math.Log2(float64(config.Ring.D)))

	t.Logf("creating public parameters")
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	A1 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	A2 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	T := fastmath.LoadNTTTransform("test", q, logD, config.Ring.BaseRing)
	p := config.P
	params := relations.KeyGenPublicParams{P1: p1, A1: A1, A2: A2, T: T, P: p}

	t.Logf("creating input")
	s := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	r := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A1, A2, s.Coeffs(), r.Coeffs(), config.P)
	k := crypto.GetAjtaiKappa(comP, comQ, config.P, config.Ring.BaseRing)
	t.Logf("generating the relation")

	rel := relations.GenerateKeyGenRelation(s, r, e, k, params, config)
	if !rel.IsValid() {
		t.Errorf("relation invalid")
	}

	t.Logf("running the protocol...")
	val := rel.CreateValidityProofDef()
	ter := rel.CreateTernaryProofDef(0, config.Ring.D)
	abp := rel.CreateApproxBoundProofDef(config.Ring.D*6, config.Ring.D*7, config.Ring.Q)
	if !fastens20.Execute(rel.S, &val, &ter, &abp, config.Ring) {
		t.Errorf("execution failed")
	}
}
