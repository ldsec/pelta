package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// RandomPoints returns a list of random points with the given mod of given size.
func RandomPoints(numPoints int, mod *big.Int, baseRing *ring.Ring) []fastmath.Coeff {
	points := make([]fastmath.Coeff, numPoints)
	for i := 0; i < numPoints; i++ {
		points[i] = fastmath.NewCoeffFromBigInt(ring.RandInt(mod), baseRing.Modulus)
	}
	return points
}

// CreateEvalMatrix constructs the evaluation matrix E where E_{i, j} = p_i^j and p_i is the ith point.
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

type Evaluator = func(p fastmath.Coeff) fastmath.Coeff

// ExtendWithPolyEval extends this linear relation system with the evaluations of the given polynomials.
func (lrbInner *LinearRelationBuilder) ExtendWithPolyEval(points []fastmath.Coeff, evaluators []Evaluator, baseRing *ring.Ring) {
	E := CreateEvalMatrix(points, baseRing)
	for i, ai := range points {
		eqnEval := NewLinearEquation(fastmath.NewIntVec(1, baseRing), 0)
		for j, evaluator := range evaluators {
			mult := evaluator(ai)
			m := fastmath.NewIntMatrixFromRows([]*fastmath.IntVec{E.RowView(i).Copy().ScaleCoeff(mult)}, baseRing)
			// TODO works ONLY for KeyGen! Fix!!
			v := lrbInner.eqns[0].rhs[j].b
			eqnEval.AppendTerm(m, v)
		}
		eqnEval.UpdateLHS()
		for j := range evaluators {
			eqnEval.AddDependency(j, j)
		}
		lrbInner.AppendEqn(eqnEval)
	}
}
