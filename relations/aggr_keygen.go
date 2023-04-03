package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"math/big"
)

func GenerateAggrKeyGenRelation(points []fastmath.Coeff, s *fastmath.Poly, r, k *fastmath.IntVec, params KeyGenPublicParams) *crypto.ImmutLinearRelation {
	lrbInner := GenerateKeyGenSystem(s, r, k, params)
	e := logging.LogExecStart("Setup.Aggregation", "extending")
	defer e.LogExecEnd()
	_, b := crypto.RLWEErrorDecomposition(fastmath.NewPoly(params.RLWEConfig.BaseRing), params.RLWEConfig)
	var evaluators []crypto.Evaluator
	evaluators = append(evaluators, func(ai fastmath.Coeff) fastmath.Coeff { return params.P1.Copy().Neg().Eval(ai) })
	for i := 0; i < b.Size(); i++ {
		iCopy := i
		evaluators = append(evaluators, func(ai fastmath.Coeff) fastmath.Coeff { return b.GetCoeff(iCopy) })
	}
	lrbInner.ExtendWithPolyEval(points, evaluators, params.RLWEConfig.BaseRing)
	return lrbInner.BuildFast(params.RLWEConfig.BaseRing)
}

func RunAggrKeyGenRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomKeyGenParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Setup.InputCreation", "working")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r, ajtaiConfig)
	k := crypto.GetAjtaiKappa(comP, comQ, ajtaiConfig)
	points := crypto.RandomPoints(4, big.NewInt(int64(rlweConfig.BaseRing.Modulus[0])), rlweConfig.BaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Setup.RelationCreation", "working")
	rel := GenerateAggrKeyGenRelation(points, s, r, k, params)
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
		WithABP(128, abpBound, fastmath.NewSlice(rlweConfig.D*6, rlweConfig.D*7)).
		WithTernarySlice(fastmath.NewSlice(0, 2*rlweConfig.D))
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
