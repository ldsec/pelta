package crypto

import "github.com/ldsec/codeBase/commitment/fastmath"

// SISProblem represents an MSIS problem, i.e., As = u for short s.
type SISProblem struct {
	A *fastmath.IntMatrix
	S *fastmath.IntVec
	U *fastmath.IntVec
}

// NewSISProblem creates a new SIS instance s.t. As = u.
func NewSISProblem(A *fastmath.IntMatrix, s *fastmath.IntVec) SISProblem {
	u := A.MulVec(s)
	return SISProblem{A, s, u}
}

// Split splits an SIS problem into multiple SIS problems over a smaller ring.
func (s SISProblem) Split(oldRing fastmath.RingParams, newRing fastmath.RingParams) []SISProblem {
	if oldRing.D <= newRing.D || oldRing.D%newRing.D != 0 {
		panic("cannot split SIS problem becasuse new ring's degree is not valid")
	}
	rebasedSIS := SISProblem{
		A: s.A.Copy().RebaseRowsLossless(newRing, 0),
		S: s.S.RebaseLossless(newRing, 0),
		U: s.U.RebaseLossless(newRing, 0),
	}
	k := len(s.S.UnderlyingPolys())
	numSplits := len(rebasedSIS.S.UnderlyingPolys()) / k
	splittedS := make([]*fastmath.IntVec, numSplits)
	for i := 0; i < numSplits; i++ {
		splitPolys := rebasedSIS.S.UnderlyingPolys()[i*k : (i+1)*k]
		splittedS[i] = fastmath.NewIntVecFromPolys(splitPolys, k*newRing.D, newRing.BaseRing)
	}
	splittedU := make([]*fastmath.IntVec, numSplits)
	j := len(s.U.UnderlyingPolys())
	for i := 0; i < numSplits; i++ {
		splitPolys := rebasedSIS.U.UnderlyingPolys()[i*j : (i+1)*j]
		splittedU[i] = fastmath.NewIntVecFromPolys(splitPolys, j*newRing.D, newRing.BaseRing)
	}
	splittedA := make([][]*fastmath.IntMatrix, numSplits)
	for i1 := 0; i1 < numSplits; i1++ {
		splittedA[i1] = make([]*fastmath.IntMatrix, numSplits)
		for i2 := 0; i2 < numSplits; i2++ {
			subsectionRows := []*fastmath.IntVec{}
			for _, r := range rebasedSIS.A.RowsView()[i1*j*newRing.D : (i1+1)*j*newRing.D] {
				relevantSubsection := r.UnderlyingPolys()[i2*k : (i2+1)*k]
				subsectionRow := fastmath.NewIntVecFromPolys(relevantSubsection, k*newRing.D, newRing.BaseRing)
				subsectionRows = append(subsectionRows, subsectionRow)
			}
			splittedA[i1][i2] = fastmath.NewIntMatrixFromRows(subsectionRows, newRing.BaseRing)
		}
	}
	splittedProblems := make([]SISProblem, 0, len(splittedA))
	for i1 := 0; i1 < numSplits; i1++ {
		for i2 := 0; i2 < numSplits; i2++ {
			splittedProblem := SISProblem{
				A: splittedA[i1][i2],
				S: splittedS[i1],
				U: splittedU[i2],
			}
			splittedProblems = append(splittedProblems, splittedProblem)
		}
	}
	return splittedProblems
}
