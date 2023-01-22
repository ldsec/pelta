package fastmath

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type PartitionedIntMatrix struct {
	rowSizes       []int
	colSizes       []int
	parts          [][]ImmutIntMatrix
	transposeCache *PartitionedIntMatrix
	intMatrixCache *IntMatrix
	baseRing       *ring.Ring
}

func NewEmptyPartitionedIntMatrix(baseRing *ring.Ring) *PartitionedIntMatrix {
	return &PartitionedIntMatrix{[]int{}, []int{}, [][]ImmutIntMatrix{}, nil, nil, baseRing}
}

func (m *PartitionedIntMatrix) Emplace(i, j int, b ImmutIntMatrix) {
	m.transposeCache = nil
	if i < len(m.parts) && m.rowSizes[i] != b.Rows() {
		panic("PartitionedIntMatrix.Emplace: row sizes do not match")
	}
	if len(m.parts) > 0 && j < len(m.parts[0]) && m.colSizes[j] != b.Cols() {
		panic("PartitionedIntMatrix.Emplace col sizes do not match")
	}
	for i >= len(m.parts) {
		upperLength := 0
		if len(m.parts) > 0 {
			upperLength = len(m.parts[0])
		}
		m.parts = append(m.parts, make([]ImmutIntMatrix, upperLength))
		m.rowSizes = append(m.rowSizes, 0)
	}
	for j >= len(m.parts[0]) {
		for i := 0; i < len(m.parts); i++ {
			m.parts[i] = append(m.parts[i], nil)
		}
		m.colSizes = append(m.colSizes, 0)
	}
	m.rowSizes[i] = b.Rows()
	m.colSizes[j] = b.Cols()
	m.parts[i][j] = b
}

func (m *PartitionedIntMatrix) PartitionHeight(i int) int {
	return m.rowSizes[i]
}

func (m *PartitionedIntMatrix) PartitionWidth(j int) int {
	return m.colSizes[j]
}

func (m *PartitionedIntMatrix) Rows() int {
	rows := 0
	for _, s := range m.rowSizes {
		rows += s
	}
	return rows
}

func (m *PartitionedIntMatrix) Cols() int {
	cols := 0
	for _, s := range m.colSizes {
		cols += s
	}
	return cols
}
func (m *PartitionedIntMatrix) String() string {
	return "PartitionedIntMatrix"
}

func (m *PartitionedIntMatrix) SizeString() string {
	return fmt.Sprintf("PartitionedIntMatrix[%d,%d]", m.Rows(), m.Cols())
}

func (m *PartitionedIntMatrix) RowsView() []*IntVec {
	return nil
}

func (m *PartitionedIntMatrix) RowView(i int) *IntVec {
	return nil
}

func (m *PartitionedIntMatrix) Hadamard(ImmutIntMatrix) ImmutIntMatrix {
	return nil
}

func (m *PartitionedIntMatrix) Transposed() ImmutIntMatrix {
	if m.transposeCache != nil {
		return logging.LogShortExecution(m.SizeString(), "transpose (cached)", func() interface{} {
			return m.transposeCache
		}).(ImmutIntMatrix)
	}
	e := logging.LogExecStart(m.SizeString(), "transpose")
	defer e.LogExecEnd()
	// Transpose the parts.
	newParts := make([][]ImmutIntMatrix, len(m.parts[0]))
	for i := 0; i < len(m.parts[0]); i++ {
		newParts[i] = make([]ImmutIntMatrix, len(m.parts))
		for j := 0; j < len(m.parts); j++ {
			if m.parts[j][i] == nil {
				continue
			}
			newParts[i][j] = m.parts[j][i].Transposed()
		}
	}
	// Transpose the size information.
	newRowSizes := make([]int, len(m.colSizes))
	newColSizes := make([]int, len(m.rowSizes))
	copy(newRowSizes, m.colSizes)
	copy(newColSizes, m.rowSizes)
	mT := &PartitionedIntMatrix{newRowSizes, newColSizes, newParts, m, nil, m.baseRing}
	m.transposeCache = mT
	return mT
}

func (m *PartitionedIntMatrix) Eq(ImmutIntMatrix) bool {
	return false
}

func (m *PartitionedIntMatrix) SubsectionCopy(int, int, int, int) *IntMatrix {
	return nil
}

func (m *PartitionedIntMatrix) GetLevel(row, col, level int) uint64 {
	return 0
}

func (m *PartitionedIntMatrix) GetCoeff(row, col int) Coeff {
	return nil
}

func (m *PartitionedIntMatrix) Max() uint64 {
	max := uint64(0)
	for _, pr := range m.parts {
		for _, p := range pr {
			if p == nil {
				continue
			} else {
				c := p.Max()
				if c > max {
					max = c
				}
			}
		}
	}
	return max
}

func (m *PartitionedIntMatrix) Copy() ImmutIntMatrix {
	c := NewEmptyPartitionedIntMatrix(m.baseRing)
	for i, pr := range m.parts {
		for j, p := range pr {
			if p == nil {
				continue
			} else {
				c.Emplace(i, j, p.Copy())
			}
		}
	}
	return c
}

func (m *PartitionedIntMatrix) MulVec(v *IntVec) *IntVec {
	if m.Cols() != v.Size() {
		panic("invalid sizes for mulvec")
	}
	e := logging.LogExecStart(fmt.Sprintf("%s.MulVec", m.SizeString()), "multiplying")
	defer e.LogExecEnd()
	out := NewIntVec(0, m.baseRing)
	for i := 0; i < len(m.parts); i++ {
		resolved := 0
		for j, p := range m.parts[i] {
			// Empty partition, skip.
			if p == nil {
				resolved += m.PartitionWidth(j)
				continue
			}
			// Apply the partition.
			vp := v.Slice(NewSlice(resolved, resolved+p.Cols()))
			out.Append(p.MulVec(vp))
			resolved += p.Cols()
		}
	}
	return out
}

func (m *PartitionedIntMatrix) AsIntMatrix() *IntMatrix {
	if m.intMatrixCache != nil {
		return m.intMatrixCache.Copy().(*IntMatrix)
	}
	var a *IntMatrix
	for i, pr := range m.parts {
		var ap *IntMatrix
		if pr[0] == nil {
			ap = NewIntMatrix(m.PartitionHeight(i), m.PartitionWidth(0), m.baseRing)
		} else {
			ap = pr[0].AsIntMatrix()
		}
		for j, p := range pr[1:] {
			if p == nil {
				ap.ExtendCols(NewIntMatrix(m.PartitionHeight(i), m.PartitionWidth(j), m.baseRing))
			} else {
				ap.ExtendCols(p)
			}
		}
		if a == nil {
			a = ap
		} else {
			a.ExtendRows(ap)
		}
	}
	m.intMatrixCache = a
	return a.Copy().(*IntMatrix)
}

func (m *PartitionedIntMatrix) RebaseRowsLossless(newRing RingParams) ImmutIntMatrix {
	e := logging.LogExecStart(m.SizeString(), "rebase")
	defer e.LogExecEnd()
	mp := NewEmptyPartitionedIntMatrix(newRing.BaseRing)
	for i, pr := range m.parts {
		for j, p := range pr {
			if p == nil {
				continue
			} else {
				mp.Emplace(i, j, p.RebaseRowsLossless(newRing))
			}
		}
	}
	return m
}

func (m *PartitionedIntMatrix) Cleanup() {
	for _, pr := range m.parts {
		for _, p := range pr {
			if p == nil {
				continue
			} else {
				p.Cleanup()
			}
		}
	}
}
