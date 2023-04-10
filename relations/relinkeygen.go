package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"math/big"
)

type RelinKeyGenPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	A  *fastmath.PolyVec
	B  *fastmath.PolyVec
	W  *fastmath.PolyVec
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
}

func GenerateRelinKeyGenRelation(s, u, er0, er1 *fastmath.Poly, r *fastmath.IntVec, params RelinKeyGenPublicParams) *crypto.ImmutLinearRelation {
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	lrb := crypto.NewLinearRelationBuilder()
	for i := 0; i < params.A.Size(); i++ {
		a := params.A.Get(i)
		b := params.B.Get(i)
		w := params.W.Get(i)
		h0i := crypto.NewLinearEquation(emptyLHS, 0).
			// aTu
			AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(a.Coeffs()), u.Coeffs()).
			// wTs
			AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(w.Coeffs()), s.Coeffs()).
			// T \sum{e0}
			AppendRLWEErrorDecompositionSum(er0, params.T, params.RLWEConfig)
		h0i.UpdateLHS()
		h1i := crypto.NewLinearEquation(emptyLHS, 0).
			// bTu
			AppendTerm(params.T.Copy().AsIntMatrix().DiagMulMat(b.Coeffs()), u.Coeffs()).
			// T \sum{e1}
			AppendRLWEErrorDecompositionSum(er1, params.T, params.RLWEConfig)
		h1i.UpdateLHS()
		// dependency on s (same s across all equations)
		if i > 0 {
			h0i.AddDependency(1, 1)
		}
		// dependency on u (h1_i depends on h0_i)
		if i == 0 {
			h1i.AddDependency(0, 0)
		} else if i == 1 {
			// first equation adds 9 independent terms
			h1i.AddDependency(0, 10)
		} else if i > 1 {
			// any new equation adds 8 independent terms
			h1i.AddDependency(0, 10+9*(i-1))
		}
		lrb.AppendEqn(h0i)
		lrb.AppendEqn(h1i)
	}
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r, params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)
	t := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r, k1, params.AjtaiConfig)
	t.AddDependency(0, 1)
	lrb.AppendEqn(t)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomRelinKeyGenParams(size int, rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) RelinKeyGenPublicParams {
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	A := fastmath.NewRandomPolyVec(size, uni, rlweConfig.BaseRing)
	B := fastmath.NewRandomPolyVec(size, uni, rlweConfig.BaseRing)
	W := fastmath.NewRandomPolyVec(size, uni, rlweConfig.BaseRing)
	T := fastmath.LoadNTTTransform("NTTTransform.test", rlweConfig.BaseRing)
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA2.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, 2*ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	params := RelinKeyGenPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		A:           A,
		B:           B,
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
	size := int((float64(logQ)+float64(logP))/float64(logP) + 0.5)
	params := getRandomRelinKeyGenParams(size, rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Setup.InputCreation", "working")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	u := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	er1 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	// k1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.RelationCreation", "working")
	rel := GenerateRelinKeyGenRelation(s, u, er0, er1, r, params)
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
		WithABP(128, abpBound, fastmath.NewSlice(rlweConfig.D*10, rlweConfig.D*11)).
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
