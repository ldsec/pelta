package fastmath

type RingElement interface {
	Eq(other RingElement) bool
}

type Vector interface {
	Get(i int) RingElement
	Set(i int, newVal RingElement) Vector
	Size() int
}

type RowMatrix interface {
	Rows() int
	Cols() int
	Get(i, j int) RingElement
	Set(i, j int, newVal RingElement)
	SetRow(i int, row Vector)
	RowView(i int) Vector
	ColCopy(i int) Vector
}

func EqVec(v1, v2 Vector) bool {
	if v1.Size() != v2.Size() {
		return false
	}
	for i := 0; i < v1.Size(); i++ {
		if !v1.Get(i).Eq(v2.Get(i)) {
			return false
		}
	}
	return true
}

func EqMatrix(m1, m2 RowMatrix) bool {
	if m1.Rows() != m2.Rows() || m1.Cols() != m2.Cols() {
		return false
	}
	for i := 0; i < m1.Rows(); i++ {
		for j := 0; j < m1.Cols(); j++ {
			if !m1.Get(i, j).Eq(m2.Get(i, j)) {
				return false
			}
		}
	}
	return true
}