package test

import (
	"math"
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/relations"
)

func getRandomKeyGenPublicParams(config relations.RelationsConfig) relations.KeyGenPublicParams {
	q := config.Ring.Q
	logD := int(math.Log2(float64(config.Ring.D)))
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	}, config.Ring.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	}, config.Ring.BaseRing)
	T := fastmath.LoadNTTTransform("test", q, logD, config.Ring.BaseRing)
	p := config.P
	params := relations.KeyGenPublicParams{P1: p1, A1: A1, A2: A2, T: T, P: p}
	return params
}

func TestKeyGen(t *testing.T) {
	config := relations.NewRelationsConfig(crypto.GetDefaultCryptoConfig())
	t.Logf("creating public parameters\n")
	params := getRandomKeyGenPublicParams(config)

	t.Logf("creating input\n")
	s := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	r := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)
	k := crypto.GetAjtaiKappa(comP, comQ, params.P, config.Ring.BaseRing)

	t.Logf("generating the relation\n")
	rel := relations.GenerateKeyGenRelation(s, r, e, k, params, config)
	if !rel.IsValid() {
		t.Errorf("relation invalid")
	}

	t.Logf("rebasing...\n")
	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebased(commitmentRing)

	t.Logf("creating protocol configuration...\n")
	protocolConfig := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel).
		WithABP(128, config.Ring.Q, fastmath.NewSlice(config.Ring.D*6, config.Ring.D*7)).
		WithTernarySlice(fastmath.NewSlice(0, config.Ring.D)).
		WithReplication(4)

	t.Logf("creating protocol parameters...\n")
	protocolParams := fastens20.GeneratePublicParameters(protocolConfig)

	t.Logf("running the protocol...\n")
	if !fastens20.Execute(rebasedRel.S, protocolParams) {
		t.Errorf("execution failed")
	}
}

func TestRelinKeyGen(t *testing.T) {
	t.Logf("reading config")
	config := relations.NewRelationsConfig(crypto.GetDefaultCryptoConfig())
	q := config.Ring.Q
	logD := int(math.Log2(float64(config.Ring.D)))

	t.Logf("creating public parameters")
	a := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	w := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	A1 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	A2 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	T := fastmath.LoadNTTTransform("test", q, logD, config.Ring.BaseRing)
	p := config.P
	l := 1
	params := relations.RelinKeyGenPublicParams{A: a, W: w, L: l, A1: A1, A2: A2, P: p, T: T}

	t.Logf("creating input")
	s := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	u := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	r := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	e0 := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	e1 := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A1, A2, s.Coeffs(), r.Coeffs(), config.P)
	k := crypto.GetAjtaiKappa(comP, comQ, config.P, config.Ring.BaseRing)

	t.Logf("generating the relation")
	rel := relations.GenerateRelinKeyGenRelation(s, u, e0, e1, r, k, params, config)
	if !rel.IsValid() {
		t.Errorf("relation invalid")
	}

	// t.Logf("running the protocol...")
	// val := rel.CreateValidityProofDef()
	// ter := rel.CreateTernaryProofDef(0, config.Ring.D)
	// abp := rel.CreateApproxBoundProofDef(config.Ring.D*6, config.Ring.D*7, config.Ring.Q)
	// if !fastens20.Execute(rel.S, &val, &ter, &abp, config.Ring) {
	// 	t.Errorf("execution failed")
	// }
}

func TestKeySwitchCollDec(t *testing.T) {
	t.Logf("reading config")
	config := relations.NewRelationsConfig(crypto.GetDefaultCryptoConfig())
	q := config.Ring.Q
	logD := int(math.Log2(float64(config.Ring.D)))

	t.Logf("creating public parameters")
	c1 := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	A1 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	A2 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	A3 := fastmath.NewIntMatrix(config.Ring.D, config.Ring.D, config.Ring.BaseRing)
	T := fastmath.LoadNTTTransform("test", q, logD, config.Ring.BaseRing)
	b := fastmath.GenerateBasis(3, config.RLWEParams.LogBeta, q, config.Ring.BaseRing)
	p := config.P
	qSmdg := big.NewInt(7)
	params := relations.KeySwitchCollDecPublicParams{C1: c1, A1: A1, A2: A2, A3: A3, T: T, B: b, P: p, QSmdg: qSmdg}

	t.Logf("creating input")
	s := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	sp := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	u := fastmath.NewRandomPoly(config.TernarySampler, config.Ring.BaseRing)
	r := fastmath.NewRandomPoly(config.UniformSampler, config.Ring.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.Ring.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A1, A2, s.Coeffs(), r.Coeffs(), config.P)
	k1 := crypto.GetAjtaiKappa(comP, comQ, config.P, config.Ring.BaseRing)
	comQ = A3.MulVec(sp.Copy().NTT().Coeffs()).Add(e.Copy().NTT().Coeffs())
	comP = comQ.Reduce(qSmdg)
	k2 := crypto.GetAjtaiKappa(comP, comQ, qSmdg, config.Ring.BaseRing)

	t.Logf("generating the relation")
	rel := relations.GenerateKeySwitchCollDecRelation(s, sp, u, r, e, k1, k2, params, config)
	if !rel.IsValid() {
		t.Errorf("relation invalid")
	}

	// t.Logf("running the protocol...")
	// val := rel.CreateValidityProofDef()
	// ter := rel.CreateTernaryProofDef(0, config.Ring.D)
	// abp := rel.CreateApproxBoundProofDef(config.Ring.D*6, config.Ring.D*7, config.Ring.Q)
	// if !fastens20.Execute(rel.S, &val, &ter, &abp, config.Ring) {
	// 	t.Errorf("execution failed")
	// }
}

func TestCollectiveBootstrapping(t *testing.T) {

}