package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type KeySwitchCollDecPublicParams struct {
	C1    *fastmath.Poly
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	T     *fastmath.IntMatrix // NTT transform
	B     *fastmath.IntVec    // ternary basis
	P     *big.Int
	QSmdg *big.Int
}

func GenerateKeySwitchCollDecRelation(s, sp, u, r, err *fastmath.Poly, k, kSmdg *fastmath.IntVec, params KeySwitchCollDecPublicParams, config RelationsConfig) crypto.LinearRelation {
	rlweParams := crypto.NewRLWEParameters(config.Ring.Q, config.Ring.D, config.Beta, config.Ring.BaseRing)
	c1T := params.C1.Coeffs().Diag().Hadamard(params.T)
	h := c1T.MulVec(s.Coeffs()).Add(err.Coeffs())
	_, t := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)
	eqn1 := crypto.NewLinearEquation(h, h.Size()).
		AppendTerm(c1T, s.Coeffs()).
		AppendVecTerm(err.Coeffs(), config.Ring.BaseRing)
	eqn2 := crypto.NewPaddedAjtaiEquation(t, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k, params.P, config.Ring.Q, config.Ring.BaseRing).
		AddDependency(0, 0)
	eqn3 := crypto.NewLinearEquation(err.Coeffs(), params.A3.Rows()).
		AppendTerm(params.A3.Copy().Hadamard(params.T), sp.Coeffs()).AddDependency(0, 0).
		AppendRLWEErrorDecompositionSum(err, params.T, rlweParams)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(eqn1).
		AppendEqn(eqn2).
		AppendEqn(eqn3)
	return lrb.Build(config.Ring.BaseRing)
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

func GeneratePubKeySwitchRelation(s, u, sp, e0, e1, r *fastmath.Poly, k1, kSmdg *fastmath.IntVec, params PubKeySwitchPublicParams, config RelationsConfig) crypto.LinearRelation {
	rlweParams := crypto.NewRLWEParameters(config.Ring.Q, config.Ring.D, config.Beta, config.Ring.BaseRing)
	c1T := params.c1.Coeffs().Diag().Hadamard(params.T)
	p0T := params.p0.Coeffs().Diag().Hadamard(params.T)
	p1T := params.p1.Coeffs().Diag().Hadamard(params.T)
	A3T := params.A3.Copy().Hadamard(params.T)
	h0 := c1T.MulVec(s.Coeffs()).Add(p0T.MulVec(u.Coeffs())).Add(e0.Coeffs())
	h1 := p1T.MulVec(u.Coeffs()).Add(e1.Copy().NTT().Coeffs())
	_, t := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.p)
	eqn1 := crypto.NewLinearEquation(h0, h0.Size()).
		AppendTerm(c1T, s.Coeffs()).
		AppendTerm(p0T, u.Coeffs()).
		AppendVecTerm(e0.Coeffs(), config.Ring.BaseRing)
	eqn2 := crypto.NewLinearEquation(h1, h1.Size()).
		AppendTerm(p1T, u.Coeffs()).AddDependency(0, 1).
		AppendRLWEErrorDecompositionSum(e1, params.T, rlweParams)
	eqn3 := crypto.NewPaddedAjtaiEquation(t, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k1, params.p, config.Ring.Q, config.Ring.BaseRing)
	eqn4 := crypto.NewLinearEquation(e0.Coeffs(), e0.Coeffs().Size()).
		AppendTerm(A3T, sp.Coeffs()).
		AppendRLWEErrorDecompositionSum(e0, params.T, rlweParams).
		AppendVecTerm(kSmdg.Copy().Scale(params.qSmdg.Uint64()).Neg(), config.Ring.BaseRing)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(eqn1).
		AppendEqn(eqn2).
		AppendEqn(eqn3).
		AppendEqn(eqn4)
	return lrb.Build(config.Ring.BaseRing)
}
