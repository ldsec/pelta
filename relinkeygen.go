package main

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type RelinKeyGenPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	A  *fastmath.PolyVec
	W  *fastmath.PolyVec
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
}

func GenerateRelinKeyGenRelation(s, u, r, er0, er1 *fastmath.Poly, k1 *fastmath.IntVec, params RelinKeyGenPublicParams) *crypto.ImmutLinearRelation {
	id := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	lrb := crypto.NewLinearRelationBuilder()
	for i := 0; i < params.A.Size(); i++ {
		a := params.A.Get(i)
		w := params.W.Get(i)
		h0i := crypto.NewLinearEquation(emptyLHS, 0).
			AppendTerm(params.T.Copy().AsIntMatrix().Neg().AsIntMatrix().DiagMulMat(a.Coeffs()), u.Coeffs()).
			AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(w.Coeffs()), s.Coeffs()).
			AppendRLWEErrorDecompositionSum(er0, params.T, params.RLWEConfig)
		h0i.UpdateLHS()
		h1i := crypto.NewLinearEquation(emptyLHS, 0).
			AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(a.Coeffs()), s.Coeffs()).
			AppendRLWEErrorDecompositionSum(er1, params.T, params.RLWEConfig)
		h1i.UpdateLHS()
		if i > 0 {
			h0iTerms := h0i.GetTerms()
			h1iTerms := h1i.GetTerms()
			for j := range h0iTerms {
				h0i.AddDependency(j, j)
			}
			for j := range h1iTerms {
				h1i.AddDependency(j, len(h0iTerms)+j)
			}
		}
		h1i.AddDependency(0, 1)
		lrb.AppendEqn(h0i)
		lrb.AppendEqn(h1i)
	}
	t := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.A1, s.Coeffs()).
		AppendTerm(params.A2, r.Coeffs()).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.AjtaiConfig.P.Uint64()), k1)
	t.UpdateLHS()
	t.AddDependency(0, 0)
	lrb.AppendEqn(t)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}
func getRandomRelinKeyGenParams(size int, rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) RelinKeyGenPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	A := fastmath.NewRandomPolyVec(size, uni, rlweConfig.BaseRing)
	W := fastmath.NewRandomPolyVec(size, uni, rlweConfig.BaseRing)
	T := fastmath.LoadNTTTransform("NTTTransform.test", rlweConfig.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	params := RelinKeyGenPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		A:           A,
		W:           W,
		A1:          A1,
		A2:          A2,
		T:           T,
	}
	return params
}
func RunRelinKeyGenRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	// (#Q+#P)/#P
	logQ := ajtaiConfig.Q.BitLen()
	logP := ajtaiConfig.P.BitLen()
	size := (logQ + logP) / logP
	params := getRandomRelinKeyGenParams(size, rlweConfig, ajtaiConfig)

	uni, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	u := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	k1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateRelinKeyGenRelation(s, u, r, er0, er1, k1, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

	e0 = logging.LogExecStart("Main", "rebasing")
	rebasedRel := rel.Rebased(rebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, rlweConfig.Q, fastmath.NewSlice(rlweConfig.D*10, rlweConfig.D*11)).
		WithTernarySlice(fastmath.NewSlice(0, 10*rlweConfig.D))
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