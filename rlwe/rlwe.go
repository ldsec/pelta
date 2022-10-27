package rlwe

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
)

// ProblemInstance represents an MRLWE problem, i.e., p0 = -p1 * s + e
type ProblemInstance struct {
	p0 rings.Polynomial
	p1 rings.Polynomial
	s  rings.Polynomial
	e  rings.Polynomial
}

func NewRLWEProblem(p1 rings.Polynomial, s rings.Polynomial, errorSampler math.PolySampler) ProblemInstance {
	e := math.NewRandomPolynomial(p1.BaseRing, errorSampler)
	p0 := p1.Copy().Neg().Mul(s).Add(e).(rings.Polynomial)
	return ProblemInstance{p0, p1, s, e}
}

// ConvertToMSIS converts an MRLWE problem to an MSIS problem in NTT domain.
func ConvertToMSIS(p ProblemInstance, nttTransform algebra.Matrix) {

}
