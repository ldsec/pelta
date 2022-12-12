package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type CollectiveBootstrappingParams struct {
	D     *fastmath.IntMatrix
	c1    *fastmath.IntVec
	a     *fastmath.IntVec
	b     *fastmath.IntVec
	A1    *fastmath.IntMatrix
	A2    *fastmath.IntMatrix
	A3    *fastmath.IntMatrix
	A4    *fastmath.IntMatrix
	p     *big.Int
	qSmdg *big.Int
	qPt   *big.Int
	T     *fastmath.IntMatrix
}

func GenerateCollectiveBootstrappingRelation(s1, s2, s3, r, e0, e1, e2 *fastmath.Poly, k1, k2, k3 *fastmath.IntVec, params CollectiveBootstrappingParams, config GlobalConfig) crypto.LinearRelation {
	rlweParams := crypto.NewRLWEParameters(config.BfvRing.Q, config.BfvRing.D, config.Beta, config.BfvRing.BaseRing)
	c1T := params.c1.Diag().Hadamard(params.T)
	negaT := params.a.Diag().Hadamard(params.T).Neg()
	M := params.A4.MulVec(s2.Copy().NTT().Coeffs()).Add(e2.Copy().NTT().Coeffs()).Add(k3.Copy().Scale(params.qPt.Uint64()).Neg())

	h0 := c1T.MulVec(s1.Coeffs()).Add(params.D.MulVec(M).Neg()).Add(e0.Coeffs())
	h1 := negaT.MulVec(s1.Coeffs()).Neg().Add(params.D.MulVec(M)).Add(e1.Copy().NTT().Coeffs())
	t, _ := crypto.GetAjtaiCommitments(params.A1, params.A2, s1.Coeffs(), r.Coeffs(), params.p)

	eqn1 := crypto.NewLinearEquation(h0, h0.Size()).
		AppendTerm(c1T, s1.Coeffs()).
		AppendTerm(params.D, M.Copy().Neg()).
		AppendVecTerm(e0.Coeffs(), config.BfvRing.BaseRing)
	eqn2 := crypto.NewLinearEquation(h1, h1.Size()).
		AppendDependentTerm(negaT.Copy(), 0).
		AppendDependentVecTerm(1, config.BfvRing.BaseRing).
		AppendRLWEErrorDecompositionSum(e1, params.T, rlweParams)
	eqn3 := crypto.NewPaddedAjtaiEquation(t, params.A1, params.A2, s1.Coeffs(), r.Coeffs(), k1, params.p, config.BfvRing.Q, config.BfvRing.BaseRing)
	eqn3.AddDependency(0, 0)
	eqn4 := crypto.NewLinearEquation(e0.Coeffs(), e0.Coeffs().Size()).
		AppendTerm(params.A3.Copy().Hadamard(params.T), s2.Coeffs()).
		AppendRLWEErrorDecompositionSum(e1, params.T, rlweParams).
		AppendVecTerm(k2.Copy().Scale(params.qSmdg.Uint64()), config.BfvRing.BaseRing)
	eqn5 := crypto.NewLinearEquation(M, M.Size()).
		AppendTerm(params.A4.Copy().Hadamard(params.T), s3.Coeffs()).
		AppendRLWEErrorDecompositionSum(e2, params.T, rlweParams).
		AppendVecTerm(k3.Copy().Scale(params.qPt.Uint64()), config.BfvRing.BaseRing)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(eqn1).
		AppendEqn(eqn2).
		AppendEqn(eqn3).
		AppendEqn(eqn4).
		AppendEqn(eqn5)
	return lrb.Build(config.BfvRing.BaseRing)
}
