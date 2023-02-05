package main

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type KeySwitchPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	C1    *fastmath.Poly
	P0p   *fastmath.Poly
	P1p   *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	QSmdg uint64
	T     *fastmath.IntMatrix // NTT transform
}

func GenerateKeySwitchRelation(s, sp, u, r, e0, e1 *fastmath.Poly, k2 *fastmath.IntVec, params KeySwitchPublicParams) *crypto.ImmutLinearRelation {
	h0 := params.C1.Copy().NTT().Mul(s.Copy().NTT()).
		Add(params.P0p.Copy().NTT().Mul(u.Copy().NTT())).
		Add(e0.Copy().NTT())
	h1 := params.P1p.Copy().NTT().Mul(u.Copy().NTT()).
		Add(e1.Copy().NTT())
	comP, comQ := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)
	eqn1 := crypto.NewLinearEquation(h0.Coeffs(), 0).
		AppendTerm(params.T.Copy().(*fastmath.IntMatrix).DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendTerm(params.T.Copy().(*fastmath.IntMatrix).DiagMulMat(params.P0p.Coeffs()), u.Coeffs()).
		AppendRLWEErrorDecompositionSum(e0, params.T, params.RLWEConfig)
	eqn2 := crypto.NewLinearEquation(h1.Coeffs(), 0).
		AppendTerm(params.T.Copy().(*fastmath.IntMatrix).DiagMulMat(params.P1p.Coeffs()), u.Coeffs()).
		AppendRLWEErrorDecompositionSum(e1, params.T, params.RLWEConfig).
		AddDependency(0, 1)
	eqn3 := crypto.NewLinearEquation(comP, 0).
		AppendTerm(params.A1, s.Coeffs()).
		AppendTerm(params.A2, r.Coeffs()).
		AppendTerm(fastmath.NewIdIntMatrix(s.N(), params.RLWEConfig.BaseRing).Scale(params.AjtaiConfig.P.Uint64()), k1).
		AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(eqn1).
		AppendEqn(eqn2).
		AppendEqn(eqn3)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomKeySwitchParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) KeySwitchPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	c1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	p0p := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	p1p := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
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
	params := KeySwitchPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          c1,
		P0p:         p0p,
		P1p:         p1p,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		QSmdg:       17,
		T:           T,
	}
	return params
}

func RunKeySwitchRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomKeySwitchParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	u := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	k2 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateKeySwitchRelation(s, sp, u, r, er0, er1, k2, params)
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
