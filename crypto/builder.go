package crypto

import (
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
	n         int
	dependent bool
}

func NewLinearEquation(lhs *fastmath.IntVec, cols int) *LinearEquation {
	return &LinearEquation{lhs, []term{}, lhs.Size(), cols, false}
}

// AppendTerm adds a new term to this equation.
func (eqn *LinearEquation) AppendTerm(A *fastmath.IntMatrix, b *fastmath.IntVec) *LinearEquation {
	if A.Rows() != eqn.m || A.Cols() != eqn.n || b.Size() != eqn.n {
		panic("cannot append term with invalid size")
	}
	eqn.rhs = append(eqn.rhs, term{A, b, false, 0})
	return eqn
}

func (eqn *LinearEquation) AppendVecTerm(b *fastmath.IntVec, baseRing *ring.Ring) *LinearEquation {
	if b.Size() != eqn.n {
		panic("cannot append term with invalid size")
	}
	if eqn.n != eqn.m {
		panic("cannot append vec term to an equation with #row != #col")
	}
	A := fastmath.NewIdIntMatrix(eqn.m, baseRing)
	eqn.rhs = append(eqn.rhs, term{A, b, false, 0})
	return eqn
}

// AddDependentTerm adds a dependent term to this equation with the associated vector being defined
// in another equation, indicated by its position in `vecIndex`.
func (eqn *LinearEquation) AppendDependentTerm(A *fastmath.IntMatrix, vecIndex int) *LinearEquation {
	eqn.rhs = append(eqn.rhs, term{A, nil, true, vecIndex})
	eqn.dependent = true
	return eqn
}

// Note: The referenced vector must have the correct size!!
func (eqn *LinearEquation) AppendDependentVecTerm(vecIndex int, baseRing *ring.Ring) *LinearEquation {
	A := fastmath.NewIdIntMatrix(eqn.m, baseRing)
	eqn.rhs = append(eqn.rhs, term{A, nil, true, vecIndex})
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

// Build constructs the linear relation of the form As = u from the appended equations.
func (lrb *LinearRelationBuilder) Build(baseRing *ring.Ring) LinearRelation {
	if len(lrb.eqns) == 0 {
		panic("cannot build a linear relation without any equations")
	}
	linRel := lrb.eqns[0].Linearize()
	for i, eqn := range lrb.eqns[1:] {
		if eqn.IsDependent() {
			prevEqnNumTerms := len(lrb.eqns[i].GetTerms())
			preB := make([]*fastmath.IntMatrix, prevEqnNumTerms)
			for _, t := range eqn.GetDependentTerms() {
				preB[t.depVecIndex] = t.A
			}
			for _, t := range eqn.GetIndependentTerms() {
				preB = append(preB, t.A)
			}
			var B *fastmath.IntMatrix
			if preB[0] == nil {
				B = fastmath.NewIntMatrix(eqn.m, eqn.n, baseRing)
			} else {
				B = preB[0]
			}
			for _, m := range preB[1:] {
				if m == nil {
					B.ExtendCols(fastmath.NewIntMatrix(eqn.m, eqn.n, baseRing))
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
