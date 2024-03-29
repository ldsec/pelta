package relations

import (
	"math/big"

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

func GenerateCollectiveBootstrappingRelation(s, sp, spp, er0, er1, er2 *fastmath.Poly, r, k2, k3 *fastmath.IntVec, params CollectiveBootstrappingPublicParams) *crypto.ImmutLinearRelation {
	id := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	DM := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().Hadamard(params.A4).AsIntMatrix().Scale(params.Delta), spp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er2, params.T.Copy().AsIntMatrix().Scale(params.Delta).AsIntMatrix(), params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.QPt*params.Delta), k3)
	DM.UpdateLHS()
	negDM := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().Hadamard(params.A4).AsIntMatrix().Scale(params.Delta).Neg(), spp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er2, params.T.Copy().AsIntMatrix().Scale(params.Delta).Neg().AsIntMatrix(), params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(params.QPt*params.Delta), k3)
	negDM.UpdateLHS()
	e0 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().Hadamard(params.A3), sp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er0, params.T, params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.QSmdg), k2)
	e0.UpdateLHS()

	h0 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendEquation(negDM).
		AppendEquation(e0)
	h0.UpdateLHS()
	h1 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.A.Coeffs()).Neg(), s.Coeffs()).
		AppendEquation(DM).
		AppendRLWEErrorDecompositionSum(er1, params.T, params.RLWEConfig)
	h1.UpdateLHS()
	// Add h1 dependencies
	h0Terms := h0.GetTerms()
	e0Terms := e0.GetTerms()
	for i := range h0Terms[:len(h0Terms)-len(e0Terms)] {
		h1.AddDependency(i, i)
	}
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r, params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)
	t := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r, k1, params.AjtaiConfig)
	t.AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(h0).
		AppendEqn(h1).
		AppendEqn(t)
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
	A2 := fastmath.PersistentIntMatrix("KeyGenA2.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, 2*ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A3 := fastmath.PersistentIntMatrix("KeyGenA3.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrixFast(rlweConfig.D, rlweConfig.D, uni, rlweConfig.BaseRing)
	}, rlweConfig.BaseRing)
	A4 := fastmath.PersistentIntMatrix("KeyGenA3.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrixFast(rlweConfig.D, rlweConfig.D, uni, rlweConfig.BaseRing)
	}, rlweConfig.BaseRing)
	delta := big.NewInt(0).Div(rlweConfig.Q, big.NewInt(int64(rlweConfig.RingParams.T)))
	params := CollectiveBootstrappingPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          c1,
		A:           a,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		A4:          A4,
		Delta:       delta.Uint64(),
		QSmdg:       uint64(1<<20) - 3,
		QPt:         rlweConfig.RingParams.T,
		T:           T,
	}
	return params
}

func RunCollectiveBootstrappingRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomCollectiveBootstrappingParams(rlweConfig, ajtaiConfig)

	uni, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Setup.InputCreation", "working")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	spp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	er2 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	// k1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	k2 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	k3 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.RelationCreation", "working")
	rel := GenerateCollectiveBootstrappingRelation(s, sp, spp, er0, er1, er2, r, k2, k3, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

	e0 = logging.LogExecStart("Setup.Rebasing", "working")
	rebasedRel := rel.Rebased(RebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.ConfigCreation", "working")
	// abp bound is set to q/2p
	abpBound := big.NewInt(0).Div(rlweConfig.Q, big.NewInt(0).Mul(ajtaiConfig.P, big.NewInt(2)))
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, abpBound, fastmath.NewSlice(rlweConfig.D*10, rlweConfig.D*12)).
		WithTernarySlice(fastmath.NewSlice(0, 10*rlweConfig.D))
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.ParamCreation", "working")
	protocolParams := fastens20.GeneratePublicParameters(protocolConfig, rebasedRel)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol execution")
	if !fastens20.Execute(rebasedRel.S, protocolParams) {
		e0.LogExecEnd()
		logging.PanicOnProduction("Main", "execution failed")
	}
	e0.LogExecEnd()
}
