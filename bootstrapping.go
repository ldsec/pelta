package main

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type CollectiveBootstrappingPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	C1    *fastmath.Poly
	A     *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	A4    *fastmath.IntMatrix
	Delta uint64
	QSmdg uint64
	QPt   uint64
	T     *fastmath.IntMatrix // NTT transform
}

func GenerateCollectiveBootstrappingRelation(s, sp, spp, r, er0, er1, er2 *fastmath.Poly, k2, k3 *fastmath.IntVec, params CollectiveBootstrappingPublicParams) *crypto.ImmutLinearRelation {
	M := params.A4.Copy().Hadamard(params.T).MulVec(spp.Copy().Coeffs()).Add(er2.Copy().NTT().Coeffs()).Add(k2.Copy().Scale(params.QPt).Neg())
	e0 := params.A3.Copy().Hadamard(params.T).MulVec(sp.Copy().Coeffs()).Add(er0.Copy().NTT().Coeffs()).Add(k3.Copy().Scale(params.QSmdg).Neg())
	h0 := params.C1.Copy().NTT().Mul(s.Copy().NTT()).Add(er0.Copy().NTT()).Add(M.Copy().Scale(params.Delta).Neg().UnderlyingPolysAsPolyNTTVec().Get(0)).Coeffs()
	h1 := params.T.Copy().AsIntMatrix().DiagMulMat(params.A.Coeffs()).Neg().MulVec(s.Coeffs()).Add(M.Copy().Scale(params.Delta)).Add(er1.Copy().NTT().Coeffs())
	comP, comQ := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)

	id := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	eqn1 := crypto.NewLinearEquation(h0, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendRLWEErrorDecompositionSum(e0.UnderlyingPolysAsPolyVec().Get(0), params.T, params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.Delta).Hadamard(params.A4).Hadamard(params.T), spp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er2, id.Copy().AsIntMatrix().Scale(-params.Delta).AsIntMatrix(), params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.Delta*params.QPt).Hadamard(params.T), k3)
	eqn2 := crypto.NewLinearEquation(h1, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.A.Coeffs()).Neg(), s.Coeffs()).AddDependency(0, 0).
		AppendRLWEErrorDecompositionSum(er1, params.T, params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(params.Delta).Hadamard(params.A4).Hadamard(params.T), spp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er2, id.Copy().AsIntMatrix().Scale(params.Delta).AsIntMatrix(), params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(params.Delta*params.QPt).Hadamard(params.T), k3)
	eqn3 := crypto.NewLinearEquation(comP, 0).
		AppendTerm(params.A1, s.Coeffs()).
		AppendTerm(params.A2, r.Coeffs()).
		AppendTerm(id.Copy().AsIntMatrix().Scale(params.AjtaiConfig.P.Uint64()), k1).
		AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder().AppendEqn(eqn1).AppendEqn(eqn2).AppendEqn(eqn3)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomCollectiveBootstrappingParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) CollectiveBootstrappingPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	c1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	a := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
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
	A4 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	params := CollectiveBootstrappingPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          c1,
		A:           a,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		A4:          A4,
		Delta:       19,
		QSmdg:       17,
		QPt:         17,
		T:           T,
	}
	return params
}

func RunCollectiveBootstrappingRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomCollectiveBootstrappingParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	spp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er2 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	k2 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing).Coeffs()
	k3 := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateCollectiveBootstrappingRelation(s, sp, spp, r, er0, er1, er2, k2, k3, params)
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
