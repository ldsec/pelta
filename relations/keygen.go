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

func GenerateKeygenRelation(s, r, err *fastmath.Poly, params KeyGenPublicParams, config GlobalConfig) crypto.LinearRelation {
	aj := crypto.NewAjtaiCommitment(params.A1, params.A2, s.Coeffs(), r.Coeffs(), config.P, config.BfvRing.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.BfvRing.Q, config.BfvRing.D, config.Beta, config.BfvRing.BaseRing)
	rlwe := crypto.NewRLWERelation(params.p1, s, err, rlweParams)
	e, b := rlwe.ErrorDecomposition()
	linRel := rlwe.ToLinearRelation(e, b, params.T)
	aj.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel
}
