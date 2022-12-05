package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type KeySwitchCollDecPublicParams struct {
	c1    *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	T     *fastmath.IntMatrix // NTT transform
	b     *fastmath.IntVec    // ternary basis
	p     *big.Int
	qSmdg *big.Int
}

func GenerateKeySwitchCollDecRelation(s, u, r, err *fastmath.Poly, k, kSmdg *fastmath.IntVec, params KeySwitchCollDecPublicParams, config GlobalConfig) crypto.LinearRelation {
	rlwe := crypto.NewRLWERelation(params.c1.Copy().Neg(), s, err, config.RLWEParams)
	aj1 := crypto.NewAjtaiCommitmentWithKappa(params.A1, params.A2, s.Coeffs(), r.Coeffs(), k, config.P)
	e, b := rlwe.ErrorDecomposition()
	linRel := rlwe.ToLinearRelation(e, b, params.T)
	aj1.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	eSum := err.SumCoeffs(0)
	aj2 := crypto.NewAjtaiCommitmentWithKappa(params.A3.Copy().Hadamard(params.T), params.T.Copy().Scale(eSum), s.Coeffs(), err.Coeffs(), kSmdg, params.qSmdg)
	aj2.EmbedIntoLinearRelation(&linRel, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel
}

type PubKeySwitchPublicParams struct {
	c1    *fastmath.Poly
	p0    *fastmath.Poly
	p1    *fastmath.Poly
	b     *fastmath.IntVec
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	p     *big.Int
	T     *fastmath.IntMatrix
	qSmdg *big.Int
}

func GeneratePubKeySwitchRelation(s, u, sp, e0, e1, r *fastmath.Poly, k1, k2 *fastmath.IntVec, params PubKeySwitchPublicParams, config GlobalConfig) crypto.LinearRelation {
	rlwe1 := crypto.NewRLWERelation(params.c1.Copy().Add(params.p0).Neg(), s.Copy().Add(u), e0, config.RLWEParams)
	rlwe2 := crypto.NewRLWERelation(params.p1, u, e1, config.RLWEParams)
	aj1 := crypto.NewAjtaiCommitment(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.p, config.BfvRing.BaseRing)
	aj2 := crypto.NewAjtaiCommitmentWithKappa(params.A3.Hadamard(params.T), params.T, s.Coeffs(), e0.Coeffs(), k2, params.qSmdg)
	aj2.ComP = e0.Coeffs()
	linRel1 := rlwe1.ToLinearRelationAuto(params.T)
	linRel2 := rlwe2.ToLinearRelationAuto(params.T)
	linRel1.AppendIndependent(linRel2)
	aj1.EmbedIntoLinearRelation(&linRel1, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	aj2.EmbedIntoLinearRelation(&linRel2, config.BfvRing.D, config.BfvRing.Q, config.BfvRing.BaseRing)
	return linRel1
}
