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
	LogD     int
	Beta     uint64
	LogBeta  int
}

func NewRLWEParameters(ens20Settings ens20.Settings) RLWEParameters {
	beta := uint64(ens20Settings.Beta())
	logBeta := int(math2.Log2(float64(beta)))
	return RLWEParameters{
		BaseRing: ens20Settings.BaseRing,
		Q:        ens20Settings.Q,
		LogD:     ens20Settings.LogD,
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
	T, _ := LoadNTTTransform(rlweProblem.Params.BaseRing, rlweProblem.Params.LogD)
	// Compute the ternary decomposition of E.
	eCoeffs := rlweProblem.E.Coeffs(rlweProblem.Params.Q)
	eDecomp, ternaryBasis := DecomposeIntoTernary(eCoeffs, rlweProblem.Params.LogBeta)
	eDecomp.Transpose()
	k := ternaryBasis.Length()
	// Compute the sub-transforms of A.
	AParts := make([]NTTTransformMatrix, k+1)
	AParts[0] = T.Extended(rlweProblem.P1.Copy().
		Neg().(rings.Polynomial).NTT())
	for i := 0; i < k; i++ {
		bi := ternaryBasis.Element(i).(rings.ZInt).Uint64()
		AParts[i+1] = T.Scaled(bi)
	}
	// Compute the sub-vectors of S.
	sParts := make([]rings.ZIntVector, k+1)
	sParts[0] = rlweProblem.S.Coeffs(rlweProblem.Params.Q)
	for i := 0; i < k; i++ {
		sParts[i+1] = rings.NewZIntVec(eDecomp.Row(i))
	}
	return SISProblemInstance{}
	// Combine.
	//A := AParts[0].AsVec()
	//for _, APart := range AParts[1:] {
	//	A = A.Appended(APart.AsVec())
	//}
	//s := sParts[0].AsVec()
	//for _, sPart := range sParts[1:] {
	//	s = s.Appended(sPart.AsVec())
	//}
	//// Compute u.
	//u := A.MulVec(s)
	//return NewSISProblem(A, rings.NewZIntVec(s), rings.NewZIntVec(u))
}
