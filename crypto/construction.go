package crypto

import "github.com/ldsec/codeBase/commitment/fastmath"

type ConstructionData struct {
	Matrices  [][]*fastmath.IntMatrix
	MatricesT [][]*fastmath.IntMatrix
}

func NewConstructionData() ConstructionData {
	return ConstructionData{Matrices: make([][]*fastmath.IntMatrix, 0)}
}

func (c ConstructionData) Insert(i, j int, matrix *fastmath.IntMatrix, matrixT *fastmath.IntMatrix) {
	// Expand as much as necessary
	for i >= len(c.Matrices) {
		c.Matrices = append(c.Matrices, []*fastmath.IntMatrix{})
		c.MatricesT = append(c.MatricesT, []*fastmath.IntMatrix{})
	}
	for k := 0; k <= i; k++ {
		for j >= len(c.Matrices[k]) {
			c.Matrices[k] = append(c.Matrices[k], nil)
			c.MatricesT[k] = append(c.MatricesT[k], nil)
		}
	}
	c.Matrices[i][j] = matrix
	c.MatricesT[i][j] = matrixT
}

func (c ConstructionData) Fill(eqns []*LinearEquation) {
	for i, eqn := range eqns {
		start := 0
		for j, term := range eqn.rhs {
			A := term.A.Copy()
			At := A.Transposed()
			if term.dependent {
				c.Insert(i, term.depVecIndex, A, At)
			} else {
				c.Insert(i, start+j, A, At)
			}
		}
		start += len(eqn.rhs)
	}
}
