package crypto

import (
	"fmt"
	"strings"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// term represents A*b
type term struct {
	A            fastmath.ImmutIntMatrix
	OriginalCols int // # of cols that are originally in this term (We maintain this to not copy A in Linearize method)
	b            *fastmath.IntVec
	dependent    bool
	depVecIndex  int
}

// LinearEquation represents a linear equation of form lhs = \sum_{i=1}^k A_i * b_i
type LinearEquation struct {
	lhs       *fastmath.IntVec
	rhs       []term
	m         int  // # rows of the equation, i.e. size of lhs.
	dependent bool // denotes whether the equation contains a dependent term.
}

// NewLinearEquation constructs an empty equation from the given lhs.
// TODO: remove the `cols` parameter
func NewLinearEquation(lhs *fastmath.IntVec, cols int) *LinearEquation {
	return &LinearEquation{lhs, []term{}, lhs.Size(), false}
}

// UpdateLHS computes the lhs of the equation and updates the `lhs` of this equation accordingly.
func (eqn *LinearEquation) UpdateLHS() *fastmath.IntVec {
	if eqn.dependent {
		panic("cannot compute dependent equation")
	}
	out := fastmath.NewIntVec(eqn.lhs.Size(), eqn.lhs.BaseRing())
	for _, t := range eqn.GetIndependentTerms() {
		out.Add(t.A.MulVec(t.b))
	}
	if !eqn.lhs.Eq(out) {
		logging.Log("LinearEquation.ComputeLHS", "computed value was different")
	}
	eqn.lhs = out
	return out
}

func (eqn *LinearEquation) GetLHS() *fastmath.IntVec {
	return eqn.lhs
}

func (eqn *LinearEquation) AppendEquation(other *LinearEquation) *LinearEquation {
	eqn.rhs = append(eqn.rhs, other.rhs...)
	return eqn
}

// AppendTerm adds a new term A*b to this equation.
func (eqn *LinearEquation) AppendTerm(A fastmath.ImmutIntMatrix, b *fastmath.IntVec) *LinearEquation {
	if A.Rows() != eqn.m || A.Cols() != b.Size() {
		panic("cannot append term with invalid size")
	}
	eqn.rhs = append(eqn.rhs, term{A, A.Cols(), b, false, 0})
	return eqn
}

// AppendVecTerm appends a new vector term Id*b to this equation.
func (eqn *LinearEquation) AppendVecTerm(b *fastmath.IntVec, baseRing *ring.Ring) *LinearEquation {
	id := fastmath.NewIdIntMatrix(b.Size(), baseRing)
	eqn.rhs = append(eqn.rhs, term{fastmath.NewCachedIntMatrix(id), id.Cols(), b, false, 0})
	return eqn
}

// AddDependentTerm adds a dependent term to this equation with the associated vector defined in another equation,
// indicated by its position in the system by `vecIndex`.
func (eqn *LinearEquation) AppendDependentTerm(A fastmath.MutIntMatrix, vecIndex int) *LinearEquation {
	eqn.rhs = append(eqn.rhs, term{A, A.Cols(), nil, true, vecIndex})
	eqn.dependent = true
	return eqn
}

// AppendDependentVecTerm adds a dependent vector term to this equation with the vecrtor defined in another equation,
// indicated by its position in the system by `vecIndex`.
// Note: The referenced vector must have the correct size (i.e., m) !!
func (eqn *LinearEquation) AppendDependentVecTerm(vecIndex int, baseRing *ring.Ring) *LinearEquation {
	id := fastmath.NewIdIntMatrix(eqn.m, baseRing)
	eqn.rhs = append(eqn.rhs, term{fastmath.NewCachedIntMatrix(id), id.Cols(), nil, true, vecIndex})
	eqn.dependent = true
	return eqn
}

// AddDependency adds a dependency to this equation.
func (eqn *LinearEquation) AddDependency(termEqnIndex, vecSystemIndex int) *LinearEquation {
	// Set the appropriate term as dependent.
	eqn.rhs[termEqnIndex].depVecIndex = vecSystemIndex
	eqn.rhs[termEqnIndex].dependent = true
	eqn.rhs[termEqnIndex].b = nil
	eqn.dependent = true
	return eqn
}

// IsDependent returns true iff this equation contains a dependent term.
func (eqn *LinearEquation) IsDependent() bool {
	return eqn.dependent
}

// GetTerms returns the list of terms in this equation.
func (eqn *LinearEquation) GetTerms() []term {
	return eqn.rhs
}

// GetDependentTerms returns the list of dependent terms in this equation.
func (eqn *LinearEquation) GetDependentTerms() []term {
	dependentTerms := []term{}
	for _, t := range eqn.GetTerms() {
		if t.dependent {
			dependentTerms = append(dependentTerms, t)
		}
	}
	return dependentTerms
}

// GetIndependentTerms returns the list of independent terms in this equation.
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
func (eqn *LinearEquation) Linearize() *LinearRelation {
	if len(eqn.rhs) == 0 {
		panic("cannot convert an empty equation into a relation")
	}
	A, ok := eqn.rhs[0].A.(*fastmath.CachedIntMatrix)
	if !ok {
		A = fastmath.NewCachedIntMatrix(eqn.rhs[0].A.AsIntMatrix())
	}
	linRel := NewLinearRelationWithLHS(A, eqn.rhs[0].b, eqn.lhs)
	for _, term := range eqn.rhs[1:] {
		linRel.ExtendPartial(term.A, term.b)
		if term.dependent {
			panic("cannot directly linearize an equation with dependent terms")
		}
	}
	return linRel
}

func (eqn *LinearEquation) SizesString() string {
	termStrs := make([]string, 0)
	for _, term := range eqn.rhs {
		var s string
		if term.dependent {
			s = fmt.Sprintf("[%d,%d][dep=%d]", term.A.Rows(), term.A.Cols(), term.depVecIndex)
		} else {
			s = fmt.Sprintf("[%d,%d][%d]", term.A.Rows(), term.A.Cols(), term.b.Size())
		}
		termStrs = append(termStrs, s)
	}
	return fmt.Sprintf("[%d] = %s", eqn.lhs.Size(), strings.Join(termStrs, " + "))
}

// String returns a string representation of this equation.
// TODO: fix
func (eqn *LinearEquation) String() string {
	return "eqn"
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

// getZeroPad constructs a zero-padding for j-th independent term of i-th equation in the given system of equations by dynamically
// selecting an appropriate column length.
func getZeroPad(i, j int, eqns []*LinearEquation, baseRing *ring.Ring) *fastmath.IntMatrix {
	// Try to find an appropriate matrix among the upper equations to determine the zero padding column size.
	curr := 0
	for _, eqn := range eqns {
		eqnIndepTerms := eqn.GetIndependentTerms()
		if j < curr+len(eqnIndepTerms) {
			// fmt.Println("found padding reference at", eqnIndex, j-curr)
			targetTerm := eqnIndepTerms[j-curr]
			return fastmath.NewIntMatrix(eqns[i].m, targetTerm.OriginalCols, baseRing)
		}
		curr += len(eqnIndepTerms)
	}
	return nil
}

// BuildFast constructs the linear relation of the form As = u from the appended equations.
func (lrb *LinearRelationBuilder) BuildFast(baseRing *ring.Ring) *ImmutLinearRelation {
	procName := "LinearRelationBuilder.BuildFast"
	if len(lrb.eqns) == 0 {
		panic("cannot build a linear relation without any equations")
	}
	e := logging.LogExecStart(procName, "building")
	defer e.LogExecEnd()
	logging.Log(procName, lrb.SizesString())
	pm := fastmath.NewEmptyPartitionedIntMatrix(baseRing)
	s := fastmath.NewIntVec(0, baseRing)
	u := fastmath.NewIntVec(0, baseRing)
	start := 0
	for i, eqn := range lrb.eqns {
		// Emplace the dependent terms.
		for _, t := range eqn.GetDependentTerms() {
			pm.Emplace(i, t.depVecIndex, t.A)
		}
		// Emplace the independent terms.
		for j, t := range eqn.GetIndependentTerms() {
			pm.Emplace(i, start+j, t.A)
			// Extend the solution s.
			s.Append(t.b)
		}
		start += len(eqn.GetIndependentTerms())
		// Extend the output.
		u.Append(eqn.lhs)
	}
	return &ImmutLinearRelation{pm, s, u}
}

// Build constructs the linear relation of the form As = u from the appended equations.
func (lrb *LinearRelationBuilder) Build(baseRing *ring.Ring) *ImmutLinearRelation {
	procName := "LinearRelationBuilder.Build"
	if len(lrb.eqns) == 0 {
		panic("cannot build a linear relation without any equations")
	}
	e := logging.LogExecStart(procName, "building")
	defer e.LogExecEnd()
	logging.Log(procName, lrb.SizesString())
	// Then construct.
	linRel := lrb.eqns[0].Linearize()
	for i, eqn := range lrb.eqns[1:] {
		// fmt.Println("LinearRelationBuilder: so far", linRel.SizesString())
		// fmt.Println("LinearRelationBuilder: built so far:", linRel.SizesString())
		if eqn.IsDependent() {
			// For a dependent equation we want to find a B, y s.t. (A || 0, B) (s, y) = (u, lhs) will yield
			// the correct linear relation.
			prevEqnNumIndepTerms := 0 // a.k.a. upper row slot amount
			for j := i; j >= 0; j-- { // i is the prev eqn index!
				prevEqnNumIndepTerms += len(lrb.eqns[j].GetIndependentTerms())
			}
			// First, emplace the dependent term matrices into the correct position.
			preB := make([]fastmath.ImmutIntMatrix, prevEqnNumIndepTerms)
			for _, t := range eqn.GetDependentTerms() {
				preB[t.depVecIndex] = t.A
			}
			// Then, append the independent term matrices.
			for _, t := range eqn.GetIndependentTerms() {
				preB = append(preB, t.A)
			}
			// fmt.Printf("lrb (%d): created B instruction of size %d\n", i+1, len(preB))
			// Iteratively build up the new row in the matrix.
			var B fastmath.MutIntMatrix
			if preB[0] == nil {
				B = getZeroPad(i+1, 0, lrb.eqns, baseRing)
			} else {
				B = preB[0].AsIntMatrix()
			}
			// fmt.Println("LinearRelationBuilder: initialized B:", B.SizeString())
			for j, m := range preB[1:] { // j is the prev independent term index (overall)!
				// fmt.Printf("lrb (%d): handling term %d\n", i+1, j+1)
				if m == nil {
					// fmt.Printf("lrb (%d): zero padding on %d %d\n", i+1, i+1, j+1)
					B.ExtendCols(getZeroPad(i+1, j+1, lrb.eqns, baseRing))
				} else {
					B.ExtendCols(m)
				}
				// fmt.Println("LinearRelationBuilder: updated B:", B.SizeString())
			}
			// Build the y, the concetanation of term vectors.
			y := eqn.GetIndependentTerms()[0].b
			for _, t := range eqn.GetIndependentTerms()[1:] {
				y.Append(t.b)
			}
			// fmt.Println("LinearRelationBuilder: extending with B:", B.SizeString())
			linRel.AppendDependent(B, y, eqn.lhs)
		} else {
			linRel.AppendIndependent(eqn.Linearize())
		}
	}
	logging.Log(procName, linRel.SizesString())
	return linRel.AsImmutable()
}

// String returns a string representation this LRB, i.e., a listing of all the appended
// equations up to now.
func (lrb *LinearRelationBuilder) String() string {
	strs := make([]string, 0, len(lrb.eqns))
	for i, eqn := range lrb.eqns {
		strs = append(strs, fmt.Sprintf("Eqn %d:\n%s", i, eqn.String()))
	}
	return strings.Join(strs, "\n")
}

// String returns a string representation this LRB, i.e., a listing of all the appended
// equations up to now.
func (lrb *LinearRelationBuilder) SizesString() string {
	strs := make([]string, 0, len(lrb.eqns))
	for i, eqn := range lrb.eqns {
		strs = append(strs, fmt.Sprintf("Eqn %d: %s", i, eqn.SizesString()))
	}
	return strings.Join(strs, "\n")
}