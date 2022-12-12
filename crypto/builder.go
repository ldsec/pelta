package crypto

import (
	"fmt"
	"strings"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type term struct {
	A           *fastmath.IntMatrix
	b           *fastmath.IntVec
	dependent   bool
	depVecIndex int
}

type LinearEquation struct {
	lhs       *fastmath.IntVec
	rhs       []term
	m         int
	dependent bool
}

func NewLinearEquation(lhs *fastmath.IntVec, cols int) *LinearEquation {
	return &LinearEquation{lhs, []term{}, lhs.Size(), false}
}

// AppendTerm adds a new term to this equation.
func (eqn *LinearEquation) AppendTerm(A *fastmath.IntMatrix, b *fastmath.IntVec) *LinearEquation {
	if A.Rows() != eqn.m || A.Cols() != b.Size() {
		panic("cannot append term with invalid size")
	}
	eqn.rhs = append(eqn.rhs, term{A, b, false, 0})
	return eqn
}

func (eqn *LinearEquation) AppendVecTerm(b *fastmath.IntVec, baseRing *ring.Ring) *LinearEquation {
	id := fastmath.NewIdIntMatrix(b.Size(), baseRing)
	eqn.rhs = append(eqn.rhs, term{id, b, false, 0})
	return eqn
}

// AddDependentTerm adds a dependent term to this equation with the associated vector being defined
// in another equation, indicated by its position in the system by `vecIndex`.
func (eqn *LinearEquation) AppendDependentTerm(A *fastmath.IntMatrix, vecIndex int) *LinearEquation {
	eqn.rhs = append(eqn.rhs, term{A, nil, true, vecIndex})
	eqn.dependent = true
	return eqn
}

// Note: The referenced vector must have the correct size!!
func (eqn *LinearEquation) AppendDependentVecTerm(vecIndex int, baseRing *ring.Ring) *LinearEquation {
	id := fastmath.NewIdIntMatrix(eqn.m, baseRing)
	eqn.rhs = append(eqn.rhs, term{id, nil, true, vecIndex})
	eqn.dependent = true
	return eqn
}

// AddDependency adds a dependency to this equation.
func (eqn *LinearEquation) AddDependency(termEqnIndex, vecSystemIndex int) *LinearEquation {
	eqn.rhs[termEqnIndex].depVecIndex = vecSystemIndex
	eqn.rhs[termEqnIndex].dependent = true
	eqn.rhs[termEqnIndex].b = nil
	eqn.dependent = true
	return eqn
}

func (eqn *LinearEquation) IsDependent() bool {
	return eqn.dependent
}

func (eqn *LinearEquation) GetTerms() []term {
	return eqn.rhs
}

func (eqn *LinearEquation) GetDependentTerms() []term {
	dependentTerms := []term{}
	for _, t := range eqn.GetTerms() {
		if t.dependent {
			dependentTerms = append(dependentTerms, t)
		}
	}
	return dependentTerms
}

func (eqn *LinearEquation) GetIndependentTerms() []term {
	independentTerms := []term{}
	for _, t := range eqn.GetTerms() {
		if !t.dependent {
			independentTerms = append(independentTerms, t)
		}
	}
	return independentTerms
}

// Linearize linearizes an independent equation.
func (eqn *LinearEquation) Linearize() LinearRelation {
	if len(eqn.rhs) == 0 {
		panic("cannot convert an empty equation into a relation")
	}
	linRel := NewLinearRelationWithLHS(eqn.rhs[0].A, eqn.rhs[0].b, eqn.lhs)
	for _, term := range eqn.rhs[1:] {
		linRel.ExtendPartial(term.A, term.b)
		if term.dependent {
			panic("cannot directly linearize an equation with dependent terms")
		}
	}
	return linRel
}

func (eqn *LinearEquation) String() string {
	strs := make([]string, 0)
	for i := 0; ; i++ {
		shouldStop := true
		rowStr := make([]string, 0)
		if eqn.lhs.Size() > i {
			rowStr = append(rowStr, fmt.Sprintf("[ %5d ]", eqn.lhs.Get(i)))
			shouldStop = false
		} else {
			rowStr = append(rowStr, "         ")
		}
		// Equal sign
		if i == eqn.lhs.Size()/2 {
			rowStr = append(rowStr, " = ")
		} else {
			rowStr = append(rowStr, "   ")
		}
		terms := make([]string, 0, len(eqn.rhs))
		for _, term := range eqn.rhs {
			termStrs := make([]string, 0, len(eqn.rhs))
			// TODO: handle empty rows
			if term.A.Rows() > i {
				termStrs = append(termStrs, term.A.RowView(i).String())
				shouldStop = false
			}
			if term.b.Size() > i {
				termStrs = append(termStrs, fmt.Sprintf("[ %5d ] ", term.b.Get(i)))
				shouldStop = false
			}
			terms = append(terms, strings.Join(termStrs, ""))
		}
		if shouldStop {
			break
		}
		termJoiner := "   "
		if eqn.lhs.Size()/2 == i {
			termJoiner = " + "
		}
		rowStr = append(rowStr, strings.Join(terms, termJoiner))
		strs = append(strs, strings.Join(rowStr, " "))
	}
	return strings.Join(strs, "\n")
}

type LinearRelationBuilder struct {
	eqns []*LinearEquation
}

// NewLinearRelationBuilder creates a new linear relation builder.
func NewLinearRelationBuilder() *LinearRelationBuilder {
	return &LinearRelationBuilder{eqns: []*LinearEquation{}}
}

// AppendEqn appends a new equation to this builder.
func (lrb *LinearRelationBuilder) AppendEqn(eqn *LinearEquation) *LinearRelationBuilder {
	lrb.eqns = append(lrb.eqns, eqn)
	return lrb
}

func getZeroPad(i, j int, eqns []*LinearEquation, baseRing *ring.Ring) *fastmath.IntMatrix {
	var targetEqn *LinearEquation
	if i > 0 {
		targetEqn = eqns[i-1]
	} else if i+1 < len(eqns) {
		targetEqn = eqns[i+1]
	} else {
		return nil
	}
	var targetMatrix *fastmath.IntMatrix
	if j < len(targetEqn.rhs) {
		targetMatrix = targetEqn.rhs[j].A
	} else {
		return nil
	}
	return fastmath.NewIntMatrix(targetMatrix.Rows(), targetMatrix.Cols(), baseRing)
}

// Build constructs the linear relation of the form As = u from the appended equations.
func (lrb *LinearRelationBuilder) Build(baseRing *ring.Ring) LinearRelation {
	if len(lrb.eqns) == 0 {
		panic("cannot build a linear relation without any equations")
	}
	linRel := lrb.eqns[0].Linearize()
	for i, eqn := range lrb.eqns[1:] {
		if eqn.IsDependent() {
			prevEqnNumIndepTerms := 0
			for j := i; j >= 0; j-- {
				prevEqnNumIndepTerms += len(lrb.eqns[j].GetIndependentTerms())
			}
			preB := make([]*fastmath.IntMatrix, prevEqnNumIndepTerms)
			for _, t := range eqn.GetDependentTerms() {
				preB[t.depVecIndex] = t.A
			}
			for _, t := range eqn.GetIndependentTerms() {
				preB = append(preB, t.A)
			}
			var B *fastmath.IntMatrix
			if preB[0] == nil {
				B = getZeroPad(i+1, 0, lrb.eqns, baseRing)
			} else {
				B = preB[0]
			}
			for j, m := range preB[1:] {
				if m == nil {
					B.ExtendCols(getZeroPad(i+1, j+1, lrb.eqns, baseRing))
				} else {
					B.ExtendCols(m)
				}
			}
			y := eqn.GetIndependentTerms()[0].b
			for _, t := range eqn.GetIndependentTerms()[1:] {
				y.Append(t.b)
			}
			linRel.AppendDependent(B, y, eqn.lhs)
		} else {
			linRel.AppendIndependent(eqn.Linearize())
		}
	}
	return linRel
}

func (lrb *LinearRelationBuilder) String() string {
	strs := make([]string, 0, len(lrb.eqns))
	for i, eqn := range lrb.eqns {
		strs = append(strs, fmt.Sprintf("Eqn %d:\n%s", i, eqn.String()))
	}
	return strings.Join(strs, "\n")
}