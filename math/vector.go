package math

import (
	"github.com/ldsec/lattigo/v2/ring"
)

type Vector struct {
	MultiArray
}

func NewVectorFromDimensions(dim int, baseRing *ring.Ring) Vector {
	return Vector{NewMultiArray([]int{dim}, baseRing)}
}

func NewVectorFromSlice(elements []*ring.Poly) Vector {
	return Vector{
		MultiArray{
			coordMap: NewCoordMap([]int{len(elements)}),
			Array:    elements,
		},
	}
}

// MapInPlace replaces every cell with the output of the given function in-place.
func (v *Vector) MapInPlace(f func(*ring.Poly, int) *ring.Poly) {
	v.MultiArray.MapInPlace(func(el *ring.Poly, coords []int) *ring.Poly {
		return f(el, coords[0])
	})
}

// ForEach calls the given function with the contents of each cell.
func (v *Vector) ForEach(f func(*ring.Poly, int)) {
	v.MultiArray.ForEach(func(el *ring.Poly, coords []int) {
		f(el, coords[0])
	})
}

// SumTo sums the polynomials in the vector and adds the result into the given `out` polynomial.
func (v *Vector) SumTo(addOp BinaryOperator, out *ring.Poly, baseRing *ring.Ring) *ring.Poly {
	v.ForEach(func(el *ring.Poly, _ int) {
		addOp.Apply(out, el, out, baseRing)
	})
	return out
}

// AddTo adds the two given vectors w.r.t. to the given addition operator and puts the result into the given
// `out` vector.
func (v *Vector) AddTo(b *Vector, addOp BinaryOperator, out *Vector, baseRing *ring.Ring) {
	// assert v.Length() == b.Length()
	for i := 0; i < v.Length(); i++ {
		addOp.Apply(v.Array[i], b.Array[i], out.Array[i], baseRing)
	}
}

// DotProductTo performs a dot product between two vectors and adds the result into the given `out` polynomial.
func (v *Vector) DotProductTo(b *Vector, mulAddOp BinaryOperator, out *ring.Poly, baseRing *ring.Ring) {
	// assert v.Length() == b.Length()
	for i := 0; i < v.Length(); i++ {
		mulAddOp.Apply(v.ElementAtIndex(i), b.ElementAtIndex(i), out, baseRing)
	}
}

// Sum sums the polynomials in the vector, returning the result.
func (v *Vector) Sum(addOp BinaryOperator, baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	v.SumTo(addOp, out, baseRing)
	return out
}

// Add adds the two given vectors w.r.t. to the given addition operator and returns the result.
func (v *Vector) Add(b *Vector, addOp BinaryOperator, baseRing *ring.Ring) Vector {
	// assert v.Length() == b.Length()
	out := NewVectorFromDimensions(v.Length(), baseRing)
	v.AddTo(b, addOp, &out, baseRing)
	return out
}

// DotProduct performs a dot product between two vectors and returns the result.
func (v *Vector) DotProduct(b *Vector, mulAddOp BinaryOperator, baseRing *ring.Ring) *ring.Poly {
	// assert v.Length() == b.Length()
	out := baseRing.NewPoly()
	v.DotProductTo(b, mulAddOp, out, baseRing)
	return out
}
