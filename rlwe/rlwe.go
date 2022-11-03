package rlwe

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

type RLWEParameters struct {
	baseRing *ring.Ring
	q        *big.Int
	logN     uint64
	beta     uint64
	logBeta  int
}

// RLWEProblemInstance represents an MRLWE problem, i.e., p0 = -p1 * s + e.
type RLWEProblemInstance struct {
	params RLWEParameters
	p0     rings.Polynomial
	p1     rings.Polynomial
	s      rings.Polynomial
	e      rings.Polynomial
}

// SISProblemInstance represents an MSIS problem, i.e., As = u for short s.
type SISProblemInstance struct {
	A algebra.Matrix
	s rings.ZIntVector
	u rings.ZIntVector
}

// NewRLWEProblem creates a new RLWE instance s.t. p0 = -p1 * s + e where e is sampled from the given error sampler.
func NewRLWEProblem(p1 rings.Polynomial, s rings.Polynomial, errorSampler math.PolySampler,
	params RLWEParameters) RLWEProblemInstance {
	e := math.NewRandomPolynomial(p1.BaseRing, errorSampler)
	p0 := p1.Copy().Neg().Mul(s).Add(e).(rings.Polynomial)
	return RLWEProblemInstance{params, p0, p1, s, e}
}

// NewSISProblem creates a new SIS instance s.t. As = u.
func NewSISProblem(A algebra.Matrix, s rings.ZIntVector, u rings.ZIntVector) SISProblemInstance {
	return SISProblemInstance{A, s, u}
}

// RLWEToSIS transforms an RLWE problem into an SIS one.
func RLWEToSIS(rlweProblem RLWEProblemInstance) SISProblemInstance {
	// Extract the NTT transform.
	T := ExtractNTTTransform(rlweProblem.params.baseRing, rlweProblem.params.logN)
	// Compute the ternary decomposition of e.
	eCoeffs := rlweProblem.e.Coeffs(rlweProblem.params.q)
	eDecomp, ternaryBasis := DecomposeIntoTernary(eCoeffs, rlweProblem.params.logBeta)
	eDecomp.Transpose()
	k := ternaryBasis.Length()
	// Compute the submatrices of A.
	AParts := make([]algebra.Matrix, k+1)
	AParts[0] = rlweProblem.p1.Copy().
		Neg().(rings.Polynomial).
		Coeffs(rlweProblem.params.q).
		Diag().
		MulMat(T)
	for i := 0; i < k; i++ {
		bi := ternaryBasis.Element(i).(rings.ZInt).Uint64()
		AParts[i+1] = T.Copy().AsMatrix().Scale(bi).AsMatrix()
	}
	// Compute the subvectors of s.
	sParts := make([]rings.ZIntVector, k+1)
	sParts[0] = rlweProblem.s.Coeffs(rlweProblem.params.q)
	for i := 0; i < k; i++ {
		sParts[i+1] = rings.NewZIntVec(eDecomp.Row(i))
	}
	// Combine.
	A := algebra.NewMatrixFromDimensions(AParts[0].Rows(), AParts[0].Cols()*(k+1))
	for i, APart := range AParts {
		APart.ForEachCol(
			func(v algebra.Vector, j int) {
				Acol := i*AParts[0].Cols() + j
				A.SetCol(Acol, v)
			})
	}
	s := algebra.NewVectorFromSize(A.Cols())
	for i, sPart := range sParts {
		sPart.ForEach(
			func(el algebra.Element, j int) {
				sCol := i*sParts[0].Length() + j
				s.SetElementAtIndex(sCol, el)
			})
	}
	// Compute u.
	u := A.MulVec(s)
	return NewSISProblem(A, rings.NewZIntVec(s), rings.NewZIntVec(u))
}
