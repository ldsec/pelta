package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

func RandomPoints(numPoints int, mod *big.Int, baseRing *ring.Ring) []fastmath.Coeff {
	points := make([]fastmath.Coeff, numPoints)
	for i := 0; i < numPoints; i++ {
		points[i] = fastmath.NewCoeffFromBigInt(ring.RandInt(mod), baseRing.Modulus)
	}
	return points
}

func CreateEvalMatrix(points []fastmath.Coeff, baseRing *ring.Ring) *fastmath.IntMatrix {
	evalMatrix := fastmath.NewIntMatrix(len(points), baseRing.N, baseRing)
	for j, p := range points {
		// let the jth row of the eval matrix to be the point's exponentiations
		evalMatrix.RowView(j).PopulateCoeffs(func(i int) fastmath.Coeff {
			out := fastmath.NewZeroCoeff(len(baseRing.Modulus))
			// compute p^i
			for lvl, qi := range baseRing.Modulus {
				out[lvl] = ring.ModExp(p[lvl], uint64(i), qi)
			}
			return out
		})
	}
	return evalMatrix
}

func (relInner *ImmutLinearRelation) ExtendWithPolyEval(numPolysPrefix int, points []fastmath.Coeff, baseRing *ring.Ring) *ImmutLinearRelation {
	eqnInner := NewLinearEquation(relInner.U, 0).AppendTerm(relInner.A, relInner.S)
	lrb2 := NewLinearRelationBuilder().AppendEqn(eqnInner)
	E := CreateEvalMatrix(points, baseRing)
	for i := 1; i < numPolysPrefix; i++ {
		E.ExtendCols(E.Copy())
	}
	AIntMatrix := relInner.A.AsIntMatrix()
	Ap := AIntMatrix.SubsectionCopy(0, numPolysPrefix*baseRing.N, 0, AIntMatrix.Cols())
	EA := E.MulMat(Ap)
	emptyLHS := fastmath.NewIntVec(EA.Rows(), baseRing)
	// TODO: precompute lhs
	eqnAggregation := NewLinearEquation(emptyLHS, 0)
	eqnAggregation.AppendTerm(EA, relInner.S)
	eqnAggregation.UpdateLHS()
	eqnAggregation.AddDependency(0, 0)
	lrb2.AppendEqn(eqnAggregation)
	return lrb2.BuildFast(baseRing)
}
