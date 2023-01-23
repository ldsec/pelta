package relations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

type KeyGenPublicParams struct {
	P1 *fastmath.Poly
	A1 *fastmath.IntMatrix
	A2 *fastmath.IntMatrix
	T  *fastmath.IntMatrix // NTT transform
	P  *big.Int
}

func GenerateKeyGenRelation(s, r, e *fastmath.Poly, k *fastmath.IntVec, params KeyGenPublicParams, config RelationsConfig) *crypto.ImmutLinearRelation {
	rlweParams := crypto.NewRLWEParameters(config.Ring.Q, config.Ring.D, config.Beta, config.Ring.BaseRing)
	_, comP := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)
	ajtaiEqn := crypto.NewPaddedAjtaiEquation(comP, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k, params.P, config.Ring.BaseRing)
	ajtaiEqn.AddDependency(0, 0)
	lrb := crypto.NewLinearRelationBuilder()
	p0 := crypto.GetRLWEP0(params.P1, s, e)
	lrb.AppendEqn(crypto.NewIndependentRLWE(p0, params.P1, s, e, params.T, rlweParams))
	lrb.AppendEqn(ajtaiEqn)
	return lrb.BuildFast(config.Ring.BaseRing)
}

// type RelinKeyGenPublicParams struct {
// 	A  *fastmath.Poly
// 	W  *fastmath.Poly
// 	L  int
// 	A1 *fastmath.IntMatrix
// 	A2 *fastmath.IntMatrix
// 	P  *big.Int
// 	T  *fastmath.IntMatrix
// }

// func GenerateRelinKeyGenRelation(s, u, e0, e1, r *fastmath.Poly, k1 *fastmath.IntVec, params RelinKeyGenPublicParams, config RelationsConfig) crypto.LinearRelation {
// 	rlweParams := crypto.NewRLWEParameters(config.Ring.Q, config.Ring.D, config.Beta, config.Ring.BaseRing)
// 	aT := params.A.Coeffs().Neg().Diag().Hadamard(params.T)
// 	wT := params.W.Coeffs().Diag().Hadamard(params.T)
// 	h0 := aT.Copy().Neg().MulVec(u.Coeffs()).Add(wT.MulVec(s.Coeffs())).Add(e0.Copy().NTT().Coeffs())
// 	h1 := aT.Copy().MulVec(s.Coeffs()).Add(e0.Copy().NTT().Coeffs())
// 	_, t := crypto.GetAjtaiCommitments(params.A1, params.A2, s.Coeffs(), r.Coeffs(), params.P)

// 	eqn1 := crypto.NewLinearEquation(h0, h0.Size()).
// 		AppendTerm(aT.Copy().Neg(), u.Coeffs()).
// 		AppendTerm(wT, s.Coeffs()).
// 		AppendRLWEErrorDecompositionSum(e0, params.T, rlweParams)
// 	eqn2 := crypto.NewLinearEquation(h1, h1.Size()).
// 		AppendDependentTerm(aT, 1).
// 		AppendRLWEErrorDecompositionSum(e1, params.T, rlweParams)
// 	eqn3 := crypto.NewPaddedAjtaiEquation(t, params.A1, params.A2, s.Coeffs(), r.Coeffs(), k1, params.P, config.Ring.BaseRing)
// 	eqn3.AddDependency(0, 1)

// 	lrb := crypto.NewLinearRelationBuilder().
// 		AppendEqn(eqn1).
// 		AppendEqn(eqn2).
// 		AppendEqn(eqn3)
// 	return lrb.Build(config.Ring.BaseRing)
// }
