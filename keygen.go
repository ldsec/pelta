package main

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type KeyGenPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	P1 *fastmath.Poly
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
}

func GenerateKeyGenRelation(s, r, e *fastmath.Poly, k *fastmath.IntVec, params KeyGenPublicParams) *crypto.ImmutLinearRelation {
	_, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.AjtaiConfig)
	ajtaiEqn := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k, params.AjtaiConfig)
	ajtaiEqn.AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder()
	p0 := crypto.RLWESample(params.P1, s, e)
	lrb.AppendEqn(crypto.NewIndependentRLWE(p0, params.P1, s, e, params.T, params.RLWEConfig))
	lrb.AppendEqn(ajtaiEqn)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomKeyGenPublicParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) KeyGenPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 0)
	p1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	T := fastmath.LoadNTTTransform("NTTTransform.test", rlweConfig.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	params := KeyGenPublicParams{P1: p1, A1: A1, A2: A2, T: T}
	return params
}

func RunKeyGenRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomKeyGenPublicParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	e := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), ajtaiConfig)
	k := crypto.GetAjtaiKappa(comP, comQ, ajtaiConfig)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateKeyGenRelation(s, r, e, k, params)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "rebasing")
	rebasedRel := rel.Rebased(rebaseRing)
	// commitmentRing := config.Ring
	// rebasedRel := rel
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, rlweConfig.Q, fastmath.NewSlice(rlweConfig.D*6, rlweConfig.D*7)).
		WithTernarySlice(fastmath.NewSlice(0, 2*rlweConfig.D)).
		WithReplication(4)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "public parameters creation")
	protocolParams := fastens20.GeneratePublicParameters(protocolConfig, rebasedRel)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol execution")
	if !fastens20.Execute(rebasedRel.S, protocolParams) {
		logging.PanicOnProduction("Main", "execution failed")
	}
	e0.LogExecEnd()
}

// type RelinKeyGenPublicParams struct {
// 	A  *fastmath.Poly
// 	W  *fastmath.Poly
// 	L  int
// 	A1 *fastmath.IntMatrix
// 	A2 *fastmath.IntMatrix
// 	P  *big.Int
// 	T  *fastmath.IntMatrix
// }

// func GenerateRelinKeyGenRelation(s, u, e0, e1, r *fastmath.Poly, k1 *fastmath.IntVec, params RelinKeyGenPublicParams, config RelationsConfig) crypto.LinearRelation {
// 	rlweParams := crypto.NewRLWEParameters(config.Ring.Q, config.Ring.D, config.Beta, config.Ring.BaseRing)
// 	aT := params.A.Coeffs().Neg().Diag().Hadamard(params.T)
// 	wT := params.W.Coeffs().Diag().Hadamard(params.T)
// 	h0 := aT.Copy().Neg().MulVec(u.Coeffs()).Add(wT.MulVec(s.Coeffs())).Add(e0.Copy().NTT().Coeffs())
// 	h1 := aT.Copy().MulVec(s.Coeffs()).Add(e0.Copy().NTT().Coeffs())
// 	_, t := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)

// 	eqn1 := crypto.NewLinearEquation(h0, h0.Size()).
// 		AppendTerm(aT.Copy().Neg(), u.Coeffs()).
// 		AppendTerm(wT, s.Coeffs()).
// 		AppendRLWEErrorDecompositionSum(e0, params.T, rlweParams)
// 	eqn2 := crypto.NewLinearEquation(h1, h1.Size()).
// 		AppendDependentTerm(aT, 1).
// 		AppendRLWEErrorDecompositionSum(e1, params.T, rlweParams)
// 	eqn3 := crypto.NewPaddedAjtaiEquation(t, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k1, params.P, config.Ring.BaseRing)
// 	eqn3.AddDependency(0, 1)

// 	lrb := crypto.NewLinearRelationBuilder().
// 		AppendEqn(eqn1).
// 		AppendEqn(eqn2).
// 		AppendEqn(eqn3)
// 	return lrb.Build(config.Ring.BaseRing)
// }
