package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type KeyGenPublicParams struct {
	p1 *fastmath.Poly
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
	b  *fastmath.IntVec    // ternary basis
	p  *big.Int
}

func GenerateKeygenRelation(s, r, err *fastmath.Poly, k *fastmath.IntVec, params KeyGenPublicParams, config GlobalConfig) crypto.LinearRelation {
	rlweParams := crypto.NewRLWEParameters(config.BfvRing.Q, config.BfvRing.D, config.Beta, config.BfvRing.BaseRing)
	_, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.p)
	ajtaiEqn := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k, params.p, config.BfvRing.Q, config.BfvRing.BaseRing)
	ajtaiEqn.AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder()
	p0 := crypto.GetRLWEP0(params.p1, s, err)
	lrb.AppendEqn(crypto.NewIndependentRLWE(p0, params.p1, s, err, params.T, rlweParams))
	lrb.AppendEqn(ajtaiEqn)
	return lrb.Build(config.BfvRing.BaseRing)
}
