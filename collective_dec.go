package main

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type CollectiveDecPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	C1    *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	QSmdg uint64
	T     *fastmath.IntMatrix // NTT transform
}

func GenerateCollectiveDecRelation(s, sp, u, r, e0 *fastmath.Poly, k2 *fastmath.IntVec, params CollectiveDecPublicParams) *crypto.ImmutLinearRelation {
	h := params.C1.Copy().NTT().Mul(s.Copy().NTT()).
		Add(params.A3.MulVec(sp.Copy().NTT().Coeffs()).UnderlyingPolysAsPolyNTTVec().Get(0)).
		Add(e0.Copy().NTT()).
		Add(k2.Copy().UnderlyingPolysAsPolyNTTVec().Get(0).Scale(params.QSmdg))
	comP, comQ := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)
	eqn1 := crypto.NewLinearEquation(h.Coeffs(), params.T.Cols()).
		AppendTerm(params.T.DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendTerm(params.T.Copy().Hadamard(params.A3), sp.Coeffs()).
		AppendRLWEErrorDecompositionSum(e0, params.T, params.RLWEConfig).
		AppendTerm(fastmath.NewIdIntMatrix(s.N(), params.RLWEConfig.BaseRing).Scale(params.QSmdg), k2)
	eqn2 := crypto.NewLinearEquation(comP, params.T.Cols()).
		AppendTerm(params.A1, s.Coeffs()).
		AppendTerm(params.A2, r.Coeffs()).
		AppendTerm(fastmath.NewIdIntMatrix(s.N(), params.RLWEConfig.BaseRing).Scale(params.AjtaiConfig.P.Uint64()), k1).
		AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(eqn1).
		AppendEqn(eqn2)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomCollectiveDecParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) CollectiveDecPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	p1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	T := fastmath.LoadNTTTransform("NTTTransform.test", rlweConfig.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A3 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	params := CollectiveDecPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          p1,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		QSmdg:       17,
		T:           T,
	}
	return params
}

func RunCollectiveDecRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomCollectiveDecParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	u := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	k2 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateCollectiveDecRelation(s, sp, u, r, er0, k2, params)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "rebasing")
	rebasedRel := rel.Rebased(rebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, rlweConfig.Q, fastmath.NewSlice(rlweConfig.D*6, rlweConfig.D*7)).
		WithTernarySlice(fastmath.NewSlice(0, 2*rlweConfig.D))
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "public parameters creation")
	protocolParams := fastens20.GeneratePublicParameters(protocolConfig, rebasedRel)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol execution")
	if !fastens20.Execute(rebasedRel.S, protocolParams) {
		e0.LogExecEnd()
		logging.PanicOnProduction("Main", "execution failed")
	}
	e0.LogExecEnd()
}
