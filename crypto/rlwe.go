package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math"
)

type RLWEParameters struct {
	BaseRing *ring.Ring
	Q        uint64
	Beta     uint64
	LogD     int
	LogBeta  int
}

func NewRLWEParameters(q uint64, logD int, beta uint64, baseRing *ring.Ring) RLWEParameters {
	logBeta := int(math.Log2(float64(beta)))
	return RLWEParameters{
		BaseRing: baseRing,
		Q:        q,
		Beta:     beta,
		LogD:     logD,
		LogBeta:  logBeta,
	}
}

// RLWEProblem represents an MRLWE problem, i.e., p0 = -p1 * s + e.
type RLWEProblem struct {
	Params RLWEParameters
	P0     *fastmath.Poly
	P1     *fastmath.Poly
	S      *fastmath.Poly
	E      *fastmath.Poly
}

// NewRLWEProblem creates a new RLWE instance s.t. p0 = -p1 * s + e where e is sampled from the given error sampler.
func NewRLWEProblem(p1 *fastmath.Poly, s *fastmath.Poly, errorSampler fastmath.PolySampler, params RLWEParameters) RLWEProblem {
	e := fastmath.NewRandomPoly(errorSampler, params.BaseRing)
	p0 := p1.Copy().NTT().Neg().Mul(s.Copy().NTT()).Add(e.Copy().NTT()).InvNTT()
	return RLWEProblem{params, p0, p1, s, e}
}

// ErrorDecomposition returns the ternary decomposition of the error.
func (r RLWEProblem) ErrorDecomposition() (*fastmath.IntMatrix, *fastmath.IntVec) {
	eCoeffs := r.E.Coeffs()
	eDecomp, ternaryBasis := fastmath.TernaryDecomposition(eCoeffs, r.Params.Beta, r.Params.LogBeta, r.Params.Q, r.Params.BaseRing)
	return eDecomp.Transposed(), ternaryBasis
}

// RLWEToSIS transforms an RLWE problem into an equivalent SIS one.
func RLWEToSIS(rlweProblem RLWEProblem) SISProblem {
	baseRing := rlweProblem.Params.BaseRing
	q := rlweProblem.Params.Q
	logD := rlweProblem.Params.LogD
	// Extract the NTT transform.
	T := fastmath.LoadNTTTransform("ntt_transform", q, logD, baseRing)
	d := T.Rows()
	// Compute the ternary decomposition of the error.
	e, b := rlweProblem.ErrorDecomposition()
	k := b.Size()
	// Compute the sub-matrices of A.
	AParts := make([]fastmath.IntMatrix, k+1)
	for i := range AParts {
		AParts[i] = *T.Copy()
	}
	negp1 := rlweProblem.P1.Copy().Neg().NTT()
	fastmath.ExtendNTTTransform(&AParts[0], negp1)
	for i := 0; i < k; i++ {
		AParts[i+1].Scale(b.Get(i))
	}
	// Compute the sub-vectors of S.
	sParts := make([]fastmath.IntVec, k+1)
	sParts[0] = *rlweProblem.S.Coeffs()
	for i := 0; i < k; i++ {
		sParts[i+1] = *e.RowView(i)
	}
	// Combine.
	A := fastmath.NewIntMatrix(d, d*(k+1), baseRing)
	A.PopulateRows(func(i int) fastmath.IntVec {
		ARowPolys := make([]fastmath.Poly, 0, k+1)
		for _, APart := range AParts {
			APartRow := APart.RowView(i)
			ARowPolys = append(ARowPolys, APartRow.UnderlyingPolys()...)
		}
		ARow := fastmath.NewIntVecFromPolys(ARowPolys, d*(k+1), baseRing)
		return *ARow
	})
	s := fastmath.NewIntVec(d*(k+1), baseRing)
	sPolys := make([]fastmath.Poly, 0, k+1)
	for _, sPart := range sParts {
		sPolys = append(sPolys, sPart.UnderlyingPolys()...)
	}
	s.SetUnderlyingPolys(sPolys)
	return NewSISProblem(A, s)
}
