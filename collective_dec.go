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
	C1    *fastmath.PolyNTT
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	QSmdg uint64
	T     *fastmath.IntMatrix // NTT transform
}

func GenerateCollectiveDecRelation(s, sp, r, er *fastmath.Poly, k1, k2 *fastmath.IntVec, params CollectiveDecPublicParams) *crypto.ImmutLinearRelation {
	id := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	e := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().Hadamard(params.A3), sp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er, params.T, params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.QSmdg), k2)
	e.UpdateLHS()
	h := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendEquation(e)
	h.UpdateLHS()
	// Create the commitment.
	t := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.A1, s.Coeffs()).
		AppendTerm(params.A2, r.Coeffs()).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.AjtaiConfig.P.Uint64()), k1)
	t.UpdateLHS()
	t.AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(h).
		AppendEqn(t)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomCollectiveDecParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) CollectiveDecPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	c1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).NTT()
	T := fastmath.LoadNTTTransform("NTTTransform.test", rlweConfig.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A3 := fastmath.PersistentIntMatrix("KeyGenA3.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrixFast(rlweConfig.D, rlweConfig.D, uni, rlweConfig.BaseRing)
	}, rlweConfig.BaseRing)
	params := CollectiveDecPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          c1,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		QSmdg:       uint64(1<<20) - 3,
		T:           T,
	}
	return params
}

func RunCollectiveDecRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomCollectiveDecParams(rlweConfig, ajtaiConfig)

	uni, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	r := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	k1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	k2 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateCollectiveDecRelation(s, sp, r, er0, k1, k2, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

	e0 = logging.LogExecStart("Main", "rebasing")
	rebasedRel := rel.Rebased(rebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, rlweConfig.Q, fastmath.NewSlice(rlweConfig.D*7, rlweConfig.D*8)).
		WithTernarySlice(fastmath.NewSlice(0, 7*rlweConfig.D))
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
