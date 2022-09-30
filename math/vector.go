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

// Sum sums the polynomials in the vector.
func (v *Vector) Sum(baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	v.ForEach(func(el *ring.Poly, _ int) {
		baseRing.Add(out, el, out)
	})
	return out
}

// Add adds the two given vectors w.r.t. to the given addition operator and returns the result.
func (v *Vector) Add(b *Vector, addOp BinaryOperator, baseRing *ring.Ring) Vector {
	// assert v.Length() == b.Length()
	out := NewVectorFromDimensions(v.Length(), baseRing)
	for i := 0; i < out.Length(); i++ {
		addOp.Apply(v.Array[i], b.Array[i], out.Array[i], baseRing)
	}
	return out
}

// AddInPlace adds the two given vectors w.r.t. to the given addition operator and adds the result
// into the given `out` vector.
func (v *Vector) AddInPlace(b *Vector, addOp BinaryOperator, out *Vector, baseRing *ring.Ring) {
	// assert v.Length() == b.Length()
	for i := 0; i < v.Length(); i++ {
		addOp.Apply(v.Array[i], b.Array[i], out.Array[i], baseRing)
	}
}

// DotProduct performs a dot product between two vectors.
func (v *Vector) DotProduct(b *Vector, mulAddOp BinaryOperator, baseRing *ring.Ring) *ring.Poly {
	out := baseRing.NewPoly()
	for i := 0; i < v.Length(); i++ {
		mulAddOp.Apply(v.ElementAtIndex(i), b.ElementAtIndex(i), out, baseRing)
	}
	return out
}

// DotProductInPlace performs a dot product between two vectors.
func (v *Vector) DotProductInPlace(b *Vector, mulAddOp BinaryOperator, out *ring.Poly, baseRing *ring.Ring) {
	for i := 0; i < v.Length(); i++ {
		mulAddOp.Apply(v.ElementAtIndex(i), b.ElementAtIndex(i), out, baseRing)
	}
}
