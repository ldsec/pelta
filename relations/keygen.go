package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type KeygenPublicParams struct {
	p1 *fastmath.Poly
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
	b  *fastmath.IntVec    // ternary basis
	p  *big.Int
}

func GenerateKeygenRelation(s, r, err *fastmath.Poly, params KeygenPublicParams, config GlobalConfig) crypto.LinearRelation {
	aj := crypto.NewAjtaiCommitment(params.A1, params.A2, s.Coeffs(), r.Coeffs(), config.P, config.BfvRing.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.BfvRing.Q, config.BfvRing.D, config.Beta, config.BfvRing.BaseRing)
	rlweRel := crypto.NewRLWERelation(params.p1, s, err, rlweParams)
	e, b := rlweRel.ErrorDecomposition()
	linRel := rlweRel.ToLinearRelation(e, b, params.T)
	aj.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel
}
