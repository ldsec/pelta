package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"math/big"
)

func GenerateAggrCollectiveDecRelation(points []fastmath.Coeff, s, sp, er *fastmath.Poly, r, k2 *fastmath.IntVec, params CollectiveDecPublicParams) *crypto.ImmutLinearRelation {
	relInner := GenerateCollectiveDecRelation(s, sp, er, r, k2, params)
	return relInner.ExtendWithPolyEval(1, points, params.RLWEConfig.BaseRing)
}

func RunAggrCollectiveDecRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomCollectiveDecParams(rlweConfig, ajtaiConfig)

	uni, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	sp := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	er0 := fastmath.NewRandomPoly(rlweConfig.ErrorSampler, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	k2 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing).Coeffs()
	points := crypto.RandomPoints(4, big.NewInt(int64(rlweConfig.BaseRing.Modulus[0])), rlweConfig.BaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateAggrCollectiveDecRelation(points, s, sp, er0, r, k2, params)
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
