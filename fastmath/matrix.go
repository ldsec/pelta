package fastmath

type ImmutIntMatrix interface {
	Rows() int
	Cols() int
	String() string
	SizeString() string
	RowsView() []*IntVec
	RowView(int) *IntVec
	MulVec(*IntVec) *IntVec
	Hadamard(ImmutIntMatrix) ImmutIntMatrix
	Transposed() ImmutIntMatrix
	Eq(ImmutIntMatrix) bool
	SubsectionCopy(int, int, int, int) *IntMatrix
	GetLevel(int, int, int) uint64
	GetCoeff(int, int) Coeff
	Max() uint64
	Copy() ImmutIntMatrix
	AsIntMatrix() *IntMatrix
	RebaseRowsLossless(RingParams) ImmutIntMatrix
	Cleanup()
}

type MutIntMatrix interface {
	ImmutIntMatrix
	ExtendRows(b ImmutIntMatrix)
	ExtendCols(b ImmutIntMatrix)
	SetForce(int, int, uint64)
	SetCoeff(int, int, Coeff)
	SetRow(int, *IntVec)
	SetCol(int, *IntVec)
	Populate(func(int, int) uint64)
	Neg() MutIntMatrix
	Scale(uint64) MutIntMatrix
	ScaleCoeff(Coeff) MutIntMatrix
}
