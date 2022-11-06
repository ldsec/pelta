package rlwe

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/ring"
	math2 "math"
	"math/big"
)

type RLWEParameters struct {
	BaseRing *ring.Ring
	Q        *big.Int
	LogD     uint64
	Beta     uint64
	LogBeta  int
}

func NewRLWEParameters(ens20Settings ens20.Settings) RLWEParameters {
	beta := uint64(ens20Settings.Beta())
	logBeta := int(math2.Log2(float64(beta)))
	return RLWEParameters{
		BaseRing: ens20Settings.BaseRing,
		Q:        ens20Settings.Q,
		LogD:     uint64(ens20Settings.LogD),
		Beta:     beta,
		LogBeta:  logBeta,
	}
}

// RLWEProblemInstance represents an MRLWE problem, i.E., P0 = -P1 * S + E.
type RLWEProblemInstance struct {
	Params RLWEParameters
	P0     rings.Polynomial
	P1     rings.Polynomial
	S      rings.Polynomial
	E      rings.Polynomial
}

// SISProblemInstance represents an MSIS problem, i.E., As = U for short S.
type SISProblemInstance struct {
	A algebra.Matrix
	S rings.ZIntVector
	U rings.ZIntVector
}

// NewRLWEProblem creates a new RLWE instance S.t. P0 = -P1 * S + E where E is sampled from the given error sampler.
func NewRLWEProblem(p1 rings.Polynomial, s rings.Polynomial, errorSampler math.PolySampler,
	params RLWEParameters) RLWEProblemInstance {
	e := math.NewRandomPolynomial(p1.BaseRing, errorSampler)
	p0 := p1.Copy().Neg().Mul(s).Add(e).(rings.Polynomial)
	return RLWEProblemInstance{params, p0, p1, s, e}
}

// NewSISProblem creates a new SIS instance S.t. As = U.
func NewSISProblem(A algebra.Matrix, s rings.ZIntVector, u rings.ZIntVector) SISProblemInstance {
	return SISProblemInstance{A, s, u}
}

// RLWEToSIS transforms an RLWE problem into an SIS one.
func RLWEToSIS(rlweProblem RLWEProblemInstance) SISProblemInstance {
	// Extract the NTT transform.
	T := GenerateNTTTransform(rlweProblem.Params.BaseRing, rlweProblem.Params.Q, rlweProblem.Params.LogD)
	// Compute the ternary decomposition of E.
	eCoeffs := rlweProblem.E.Coeffs(rlweProblem.Params.Q)
	eDecomp, ternaryBasis := DecomposeIntoTernary(eCoeffs, rlweProblem.Params.LogBeta)
	eDecomp.Transpose()
	k := ternaryBasis.Length()
	// Compute the submatrices of A.
	AParts := make([]algebra.Matrix, k+1)
	AParts[0] = rlweProblem.P1.Copy().
		Neg().(rings.Polynomial).
		Coeffs(rlweProblem.Params.Q).
		Diag().
		MulMat(T)
	for i := 0; i < k; i++ {
		bi := ternaryBasis.Element(i).(rings.ZInt).Uint64()
		AParts[i+1] = T.Copy().AsMatrix().Scale(bi).AsMatrix()
	}
	// Compute the subvectors of S.
	sParts := make([]rings.ZIntVector, k+1)
	sParts[0] = rlweProblem.S.Coeffs(rlweProblem.Params.Q)
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
	// Compute U.
	u := A.MulVec(s)
	return NewSISProblem(A, rings.NewZIntVec(s), rings.NewZIntVec(u))
}
