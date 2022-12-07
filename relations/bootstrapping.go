package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type CollectiveBootstrappingParams struct {
	D     *fastmath.IntVec
	c1    *fastmath.IntVec
	a     *fastmath.IntVec
	b     *fastmath.IntVec
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	p     *big.Int
	qSmdg *big.Int
	qPt   *big.Int
	T     *fastmath.IntMatrix
}

func GenerateCollectiveBootstrappingRelation(s1, s2, s3, e0, e1, e2, r, k1, k2, k3 *fastmath.IntVec, params CollectiveBootstrappingParams, config GlobalConfig) crypto.LinearRelation {
	lrb := crypto.NewLinearRelationBuilder()
	h0 := params.c1.Diag().Hadamard(params.T).MulVec(s1).
		Add(params.D.Copy().Neg()).
		Add(e0)
	lrb.AppendEqn(
		crypto.NewLinearEquation(h0, h0.Size()).
			AppendTerm(params.c1.Diag().Hadamard(params.T), s1).
			AppendVecTerm(params.D, config.BfvRing.BaseRing).
			AppendVecTerm(e0, config.BfvRing.BaseRing))
	// h1 := params.a.Copy().Neg().Diag().Hadamard(params.T).MulVec(s1).Add(params.D)
	// lrb.AppendEqn()

	return lrb.Build(config.BfvRing.BaseRing)
}
