package relations

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

type AggregationPublicParams struct {
	crypto.AjtaiConfig
	crypto.RLWEConfig
	Points []uint64
	A1     *fastmath.IntMatrix
	A2     *fastmath.IntMatrix
}

func GenerateAggregationRelation(s *fastmath.Poly, r *fastmath.IntVec, params AggregationPublicParams) *crypto.ImmutLinearRelation {
	// evaluation equations
	eqns := make([]*crypto.LinearEquation, len(params.Points))
	emptyLHS := fastmath.NewIntVec(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
	for j, p := range params.Points {
		// create the evaluation matrix
		evalMatrix := fastmath.NewIdIntMatrix(params.RLWEConfig.D, params.RLWEConfig.BaseRing)
		for i := 0; i < params.RLWEConfig.D; i++ {
			// compute p to the ith power
			pi := big.NewInt(0).Exp(big.NewInt(int64(p)), big.NewInt(int64(i)), params.RLWEConfig.Q)
			// coeffs
			piCoeffs := fastmath.NewCoeffFromBigInt(pi, params.RLWEConfig.BaseRing.Modulus)
			evalMatrix.RowView(i).ScaleCoeff(piCoeffs)
		}
		// add the evaluation equation
		eqns[j] = crypto.NewLinearEquation(emptyLHS, 0).AppendTerm(evalMatrix, s.Coeffs())
		eqns[j].UpdateLHS()
		if j > 0 {
			eqns[j].AddDependency(0, 0)
		}
	}
	// commitment equation
	comQ, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r, params.AjtaiConfig)
	k1 := crypto.GetAjtaiKappa(comP, comQ, params.AjtaiConfig)
	t := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r, k1, params.AjtaiConfig)
	t.AddDependency(0, 0)
	// create the relation
	lrb := crypto.NewLinearRelationBuilder()
	for _, eqn := range eqns {
		lrb.AppendEqn(eqn)
	}
	lrb.AppendEqn(t)
	return lrb.BuildFast(params.RLWEConfig.BaseRing)
}

func getRandomAggregationParams(numPoints int, rlweConfig crypto.RLWEConfig, ajtaiConfig crypto.AjtaiConfig) AggregationPublicParams {
	A1 := fastmath.PersistentIntMatrix("KeyGenA1.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	A2 := fastmath.PersistentIntMatrix("KeyGenA2.test", func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrix(ajtaiConfig.D, 2*ajtaiConfig.D, ajtaiConfig.P, ajtaiConfig.BaseRing)
	}, ajtaiConfig.BaseRing)
	// Get the smallest modulus q_i. The points must be sampled from Z_{q_i}
	qi := rlweConfig.RingParams.BaseRing.Modulus[0]
	points := make([]uint64, numPoints)
	for i := 0; i < numPoints; i++ {
		points[i] = ring.RandInt(big.NewInt(int64(qi))).Uint64()
	}
	return AggregationPublicParams{
		AjtaiConfig: ajtaiConfig,
		RLWEConfig:  rlweConfig,
		A1:          A1,
		A2:          A2,
		Points:      points,
	}
}

func RunAggregationRelation() {
	rlweConfig := GetDefaultRLWEConfig()
	ajtaiConfig := GetDefaultAjtaiConfig()
	params := getRandomAggregationParams(1, rlweConfig, ajtaiConfig)

	_, ter, _ := fastmath.GetSamplers(rlweConfig.RingParams, 1)
	e0 := logging.LogExecStart("Main", "input creation")
	s := fastmath.NewRandomPoly(ter, rlweConfig.BaseRing)
	r := fastmath.NewRandomIntVecFast(params.A2.Cols(), ter, rlweConfig.BaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "relation creation")
	rel := GenerateAggregationRelation(s, r, params)
	e0.LogExecEnd()
	if !rel.IsValid() {
		panic("invalid relation")
	}

	e0 = logging.LogExecStart("Main", "rebasing")
	rebasedRel := rel.Rebased(rebaseRing)
	e0.LogExecEnd()

	e0 = logging.LogExecStart("Main", "protocol config creation")
	protocolConfig := GetDefaultProtocolConfig(rebasedRel.A.Rows(), rebasedRel.A.Cols()).
		WithABP(128, rlweConfig.Q, fastmath.NewSlice(rlweConfig.D*10, rlweConfig.D*12)).
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
