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
	rlwe := crypto.NewRLWERelation(params.p1, s, err, rlweParams)
	e, b := rlwe.ErrorDecomposition()
	linRel := rlwe.ToLinearRelation(e, b, params.T)
	aj.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel
}

type KeyswitchPublicParams struct {
	c1    *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	T     *fastmath.IntMatrix // NTT transform
	b     *fastmath.IntVec    // ternary basis
	p     *big.Int
	qSmdg *big.Int
}

func GenerateKeyswitchRelation(s, u, err, r *fastmath.Poly, k, kSmdg uint64, params KeyswitchPublicParams, config GlobalConfig) crypto.LinearRelation {
	aj1 := crypto.NewAjtaiCommitment(params.A1, params.A2, s.Coeffs(), r.Coeffs(), config.P, config.BfvRing.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.BfvRing.Q, config.BfvRing.D, config.Beta, config.BfvRing.BaseRing)
	rlwe := crypto.NewRLWERelation(params.c1.Copy().Neg(), s, err, rlweParams)
	e, b := rlwe.ErrorDecomposition()
	linRel := rlwe.ToLinearRelation(e, b, params.T)
	aj1.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	eSum := err.SumCoeffs(0)
	aj2 := crypto.NewAjtaiCommitment(params.A3.Copy().Hadamard(params.T), params.T.Copy().Scale(eSum), s.Coeffs(), err.Coeffs(), params.qSmdg, config.BfvRing.BaseRing)
	aj2.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel
}