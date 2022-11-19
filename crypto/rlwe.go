package crypto

import (
	"math"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type RLWEParameters struct {
	BaseRing *ring.Ring
	Q        uint64
	Beta     uint64
	LogD     int
	LogBeta  int
}

func NewRLWEParameters(q uint64, d int, beta uint64, baseRing *ring.Ring) RLWEParameters {
	logBeta := int(math.Log2(float64(beta)))
	logD := int(math.Log2(float64(d)))
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
	// Compute the ternary decomposition of the error.
	e, b := rlweProblem.ErrorDecomposition()
	k := b.Size()
	// Compute the sub-matrices of A.
	aParts := make([]*fastmath.IntMatrix, k+1)
	for i := range aParts {
		aParts[i] = T.Copy()
	}
	negP1 := rlweProblem.P1.Copy().Neg().NTT()
	fastmath.ExtendNTTTransform(aParts[0], negP1)
	for i := 0; i < k; i++ {
		aParts[i+1].Scale(b.Get(i))
	}
	// Compute the sub-vectors of S.
	sParts := make([]*fastmath.IntVec, k+1)
	sParts[0] = rlweProblem.S.Coeffs()
	for i := 0; i < k; i++ {
		sParts[i+1] = e.RowView(i)
	}
	// Combine.
	A := aParts[0]
	for _, aPart := range aParts[1:] {
		A.ExtendCols(aPart)
	}
	s := sParts[0]
	for _, sPart := range sParts[1:] {
		s.Append(sPart)
	}
	return NewSISProblem(A, s)
}
