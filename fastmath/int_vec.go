package fastmath

import (
	"fmt"
	"math/big"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

type IntVec struct {
	size     int
	polys    []*Poly
	baseRing *ring.Ring
	// set to the version of this vector before rebasing
	unrebasedRef *IntVec
}

func NewIntVec(size int, baseRing *ring.Ring) *IntVec {
	numPolys := int(size/baseRing.N) + 1
	if size%baseRing.N == 0 {
		numPolys -= 1
	}
	polys := make([]*Poly, numPolys)
	for i := 0; i < len(polys); i++ {
		polys[i] = NewPoly(baseRing)
	}
	return &IntVec{size, polys, baseRing, nil}
}

func NewIntVecFromSlice(slice []uint64, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(len(slice), baseRing)
	for i := 0; i < len(slice); i++ {
		v.SetForce(i, slice[i])
	}
	return v
}
func NewIntVecFromCoeffSlice(slice []Coeff, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(len(slice), baseRing)
	for i := 0; i < len(slice); i++ {
		v.SetCoeff(i, slice[i])
	}
	return v
}

func NewIntVecFromPolys(polys []*Poly, size int, baseRing *ring.Ring) *IntVec {
	return &IntVec{size, polys, baseRing, nil}
}

// Size returns the size of this vector.
func (v *IntVec) Size() int {
	return v.size
}

// Populate is used to initialize the elements of this vector.
func (v *IntVec) Populate(f func(int) uint64) {
	for i := 0; i < v.size; i++ {
		val := f(i)
		v.SetForce(i, val)
	}
}

// Populate is used to initialize the elements of this vector.
func (v *IntVec) PopulateCoeffs(f func(int) Coeff) {
	for i := 0; i < v.size; i++ {
		val := f(i)
		v.SetCoeff(i, val)
	}
}

// GetCoeff returns the element at the given index
func (v *IntVec) GetCoeff(index int) Coeff {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	return v.polys[polyIndex].GetCoeff(coeffIndex)
}

// GetLevel returns the element at the given index and level
func (v *IntVec) GetLevel(index, level int) uint64 {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	return v.polys[polyIndex].GetLevel(coeffIndex, level)
}

// UnderlyingPolys returns the polynomials that are being used to represent this vector.
func (v *IntVec) UnderlyingPolys() []*Poly {
	return v.polys
}

func (v *IntVec) UnderlyingPolysAsIntVecs() []*IntVec {
	vecs := make([]*IntVec, len(v.polys))
	for i, p := range v.polys {
		vecs[i] = NewIntVecFromPolys([]*Poly{p}, p.N(), p.baseRing)
	}
	if v.Size()%v.baseRing.N != 0 {
		vecs[len(vecs)-1].size = v.Size() % v.baseRing.N
	}
	return vecs
}

func (v *IntVec) UnderlyingPolysAsPolyVec() *PolyVec {
	pV := NewPolyVec(len(v.polys), v.BaseRing())
	for i, p := range v.polys {
		pV.Set(i, p)
	}
	return pV
}

func (v *IntVec) UnderlyingPolysAsPolyNTTVec() *PolyNTTVec {
	pV := NewPolyVec(len(v.polys), v.BaseRing()).NTT()
	for i, p := range v.polys {
		pV.Set(i, ForceNTT(p))
	}
	return pV
}

// SetUnderlyingPolys can be used to update the polynomials that are used to represent this vector.
func (v *IntVec) SetUnderlyingPolys(polys []*Poly) {
	v.polys = polys
}

// SetForce updates the given element of this vector.
func (v *IntVec) SetForce(index int, newValue uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].SetForce(coeffIndex, newValue)
}

func (v *IntVec) SetCoeff(index int, newValues []uint64) {
	polyIndex := index / v.baseRing.N
	coeffIndex := index % v.baseRing.N
	v.polys[polyIndex].SetCoeff(coeffIndex, newValues)
}

// Copy copies the vector.
func (v *IntVec) Copy() *IntVec {
	polys := make([]*Poly, len(v.polys))
	for i, p := range v.polys {
		polys[i] = p.Copy()
	}
	unrebased := v.unrebasedRef
	// if unrebased != nil {
	// 	unrebased = unrebased.Copy()
	// }
	return &IntVec{v.size, polys, v.baseRing, unrebased}
}

// String returns the string representation of this integer vector.
func (v *IntVec) String() string {
	truncatedSize := 10
	if v.Size() < truncatedSize {
		truncatedSize = v.Size()
	}
	s := fmt.Sprintf("IntVec[%d]{", v.Size())
	elemStrs := make([]string, 0, v.Size())
	for i := 0; i < v.Size(); i++ {
		elemStrs = append(elemStrs, fmt.Sprintf("%v", v.GetCoeff(i)))
	}
	return s + strings.Join(elemStrs[:truncatedSize], ",") + "}"
}

// Append appends the contents of the given vector into this one.
func (v *IntVec) Append(r *IntVec) *IntVec {
	// Update the size.
	newSize := v.size + r.size
	// Optimization.
	if v.Size()%v.baseRing.N == 0 {
		v.size = newSize
		v.polys = append(v.polys, r.polys...)
		return v
	}
	// Extend the number of underlying polynomials.
	numPolys := int(newSize/v.baseRing.N) + 1
	if newSize%v.baseRing.N == 0 {
		numPolys -= 1
	}
	oldSize := v.size
	v.size = newSize
	for i := len(v.polys); i < numPolys; i++ {
		v.polys = append(v.polys, NewPoly(v.baseRing))
	}
	// Move the elements in.
	for i := 0; i < r.Size(); i++ {
		v.SetCoeff(oldSize+i, r.GetCoeff(i))
	}
	return v
}

// SliceCopy returns a (copied) slice of this vector.
func (v *IntVec) SliceWithPolys(startPoly, endPoly, size int) *IntVec {
	return NewIntVecFromPolys(v.polys[startPoly:endPoly], size, v.baseRing)
}

func (v *IntVec) Cleanup() {
	for _, p := range v.polys {
		p.Cleanup()
	}
}

// RebaseLossless rebases every underlying polynomial, splitting them when necessary so that
// no coefficient will be discarded.
func (v *IntVec) RebaseLossless(newRing RingParams) *IntVec {
	if v.baseRing.N <= newRing.D || v.baseRing.N%newRing.D != 0 {
		panic(fmt.Sprintf("cannot rebase lossless %d to %d", v.baseRing.N, newRing.D))
	}
	// Each underlying polynomial will be represented by this many rebased polynomials.
	splitsPerPoly := v.baseRing.N / newRing.D
	newPolys := make([]*Poly, 0, splitsPerPoly*len(v.polys))
	for _, p := range v.polys {
		newPolys = append(newPolys, p.SplitCoeffs(newRing.BaseRing)...)
	}
	vp := NewIntVecFromPolys(newPolys, v.Size(), newRing.BaseRing)
	vp.unrebasedRef = v
	return vp
}

func (v *IntVec) All(pred func(el uint64) bool) bool {
	for i := 0; i < v.Size(); i++ {
		if !pred(v.GetLevel(i, 0)) {
			return false
		}
	}
	return true
}

// BaseRing returns the polynomial ring over which this integer vector is defined.
func (v *IntVec) BaseRing() *ring.Ring {
	return v.baseRing
}

func (v *IntVec) GetAllLevel(level int) []uint64 {
	// a := make([]uint64, len(v.polys[0].ref.Coeffs[level]))
	// for i := 0; i < len(a); i++ {
	// 	a[i] = v.polys[0].ref.Coeffs[level][i]
	// }
	a := v.polys[0].ref.Coeffs[level]
	for _, b := range v.polys[1:] {
		a = append(a, b.ref.Coeffs[level]...)
	}
	return a[:v.Size()]
}

func (v *IntVec) SetAllLevel(level int, vals []uint64) {
	start := 0
	for _, p := range v.polys {
		end := start + p.N()
		if end > len(vals) {
			end = len(vals)
			for i := start; i < end; i++ {
				p.SetLevel(i, level, vals[i])
			}
			return
		}
		p.SetAllLevel(level, vals[start:end])
		start += p.N()
	}
}

func (v *IntVec) IsUnset() bool {
	for _, p := range v.polys {
		if !p.IsUnset() {
			return false
		}
	}
	return true
}

// Reduce reduces the elements of this vector by the given mod.
func (v *IntVec) Reduce(mod *big.Int) *IntVec {
	for _, p := range v.polys {
		for i := 0; i < p.N(); i++ {
			for lvl := 0; lvl < len(p.ref.Coeffs); lvl++ {
				reducedVal := big.NewInt(0).Mod(big.NewInt(int64(p.GetLevel(i, lvl))), mod)
				p.SetLevel(i, lvl, reducedVal.Uint64())
			}
		}
	}
	return v
}

// Scale scales the vector by the given scalar factor amount.
func (v *IntVec) Scale(factor uint64) *IntVec {
	for _, p := range v.polys {
		p.Scale(factor)
	}
	return v
}

// ScaleCoeff scales the vector by the given coefficient.
func (v *IntVec) ScaleCoeff(factors Coeff) *IntVec {
	for _, p := range v.polys {
		p.ScaleCoeff(factors)
	}
	return v
}

// Neg negates the polynomial.
func (v *IntVec) Neg() *IntVec {
	for _, p := range v.polys {
		p.Neg()
	}
	return v
}

// Diag returns a diagonalization of this vector.
func (v *IntVec) Diag() *IntMatrix {
	m := NewIntMatrix(v.Size(), v.Size(), v.BaseRing())
	for i := 0; i < v.Size(); i++ {
		m.SetCoeff(i, i, v.GetCoeff(i))
	}
	return m
}

// Max returns the largest element of the vector.
func (v *IntVec) Max(q *big.Int) *big.Int {
	max := v.polys[0].Max(q)
	for _, p := range v.polys[1:] {
		pMax := p.Max(q)
		if pMax.Cmp(max) > 0 {
			max = pMax
		}
	}
	return max
}

// Add adds two integer vectors.
func (v *IntVec) Add(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Add sizes do not match")
	}
	numPolys := len(v.polys)
	if len(r.polys) < numPolys {
		numPolys = len(r.polys)
	}
	for i := 0; i < numPolys; i++ {
		v.polys[i].Add(r.polys[i])
	}
	return v
}

func (v *IntVec) Sum() Coeff {
	sum := NewPoly(v.baseRing)
	for _, p := range v.polys {
		sum.Add(p)
	}
	return sum.SumCoeffs()
}

// Dot returns the dot product of the given two vectors.
func (v *IntVec) Dot(r *IntVec) Coeff {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	// if the vectors are both rebased, do the dot product with the unrebased version
	if v.unrebasedRef != nil && r.unrebasedRef != nil {
		return v.unrebasedRef.Dot(r.unrebasedRef)
	}
	numPolys := len(v.polys)
	if len(r.polys) < numPolys {
		numPolys = len(r.polys)
	}
	preSum := NewPoly(v.baseRing)
	for i := 0; i < numPolys; i++ {
		a := v.polys[i]
		b := r.polys[i]
		if a.IsUnset() || b.IsUnset() {
			continue
		}
		a.MulCoeffsAndAdd(b, preSum)
	}
	return preSum.SumCoeffs()
}

// Hadamard performs coefficient-wise multiplication.
func (v *IntVec) Hadamard(r *IntVec) *IntVec {
	if v.size != r.size {
		panic("IntVec.Dot sizes do not match")
	}
	numPolys := len(v.polys)
	if len(r.polys) < numPolys {
		numPolys = len(r.polys)
	}
	for i := 0; i < numPolys; i++ {
		a := v.polys[i]
		b := r.polys[i]
		a.MulCoeffs(b)
	}
	return v
}

// Eq checks the equality between two integer vectors.
func (v *IntVec) Eq(r *IntVec) bool {
	if v.Size() != r.Size() {
		return false
	}
	// The underlying polynomial # might change even when the sizes are equal.
	// Take the minimum.
	numPolys := len(v.polys)
	if len(r.polys) < numPolys {
		numPolys = len(r.polys)
	}
	for i := 0; i < numPolys; i++ {
		if !v.polys[i].Eq(r.polys[i]) {
			return false
		}
	}
	return true
}
