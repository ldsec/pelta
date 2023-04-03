package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"math/big"
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

func GenerateKeySwitchRelation(s, sp, u, er0, er1 *fastmath.Poly, r, k2 *fastmath.IntVec, params KeySwitchPublicParams) *crypto.ImmutLinearRelation {
	id := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	e0 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.A3.Copy().Hadamard(params.T), sp.Coeffs()).
		AppendRLWEErrorDecompositionSum(er0, params.T, params.RLWEConfig).
		AppendTerm(id.Copy().AsIntMatrix().Scale(-params.QSmdg), k2)
	e0.UpdateLHS()
	h0 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.C1.Coeffs()), s.Coeffs()).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.P0p.Coeffs()), u.Coeffs())
	h0.UpdateLHS()
	h0.AppendEquation(e0)
	h1 := crypto.NewLinearEquation(emptyLHS, 0).
		AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(params.P1p.Coeffs()), u.Coeffs()).
		AppendRLWEErrorDecompositionSum(er1, params.T, params.RLWEConfig)
	h1.UpdateLHS()
	h1.AddDependency(0, 1)
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

func getRandomKeySwitchParams(rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) KeySwitchPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	c1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	p0p := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	p1p := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
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
	params := KeySwitchPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		C1:          c1,
		P0p:         p0p,
		P1p:         p1p,
		A1:          A1,
		A2:          A2,
		A3:          A3,
		QSmdg:       uint64(1<<20) - 3,
		T:           T,
	}
	return params
}

func RunKeySwitchRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomKeySwitchParams(rlweConfig, ajtaiConfig)

	uni, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Setup.InputCreation", "working")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	u := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	// k1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	k2 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.RelationCreation", "working")
	rel := GenerateKeySwitchRelation(s, sp, u, er0, er1, r, k2, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

	e0 = logging.LogExecStart("Setup.Rebasing", "working")
	rebasedRel := rel.Rebased(rebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.ConfigCreation", "working")
	// abp bound is set to q/2p
	abpBound := big.NewInt(0).Div(rlweConfig.Q, big.NewInt(0).Mul(ajtaiConfig.P, big.NewInt(2)))
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, abpBound, fastmath.NewSlice(rlweConfig.D*11, rlweConfig.D*13)).
		WithTernarySlice(fastmath.NewSlice(0, 11*rlweConfig.D))
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
