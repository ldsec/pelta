package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"math/big"
)

func GenerateAggrKeyGenRelation(points []fastmath.Coeff, s *fastmath.Poly, r, k *fastmath.IntVec, params KeyGenPublicParams) *crypto.ImmutLinearRelation {
	relInner := GenerateKeyGenRelation(s, r, k, params)
	return relInner.ExtendWithPolyEval(1, points, params.RLWEConfig.BaseRing)
}

func RunAggrKeyGenRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomKeyGenParams(rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	points := crypto.RandomPoints(4, big.NewInt(int64(rlweConfig.BaseRing.Modulus[0])), rlweConfig.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r, ajtaiConfig)
	k := crypto.GetAjtaiKappa(comP, comQ, ajtaiConfig)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateAggrKeyGenRelation(points, s, r, k, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

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
