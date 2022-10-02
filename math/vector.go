package math

import (
	"github.com/ldsec/lattigo/v2/ring"
)

type Vector struct {
	*MultiArray
}

// NewVectorFromSize constructs a new empty vector with the given size.
func NewVectorFromSize(dim int, baseRing *ring.Ring) Vector {
	return Vector{NewMultiArray([]int{dim}, baseRing)}
}

// NewVectorFromSlice constructs a new vector from the given slice.
// Warning: Does not copy the underlying elements.
func NewVectorFromSlice(elements []RingElement, baseRing *ring.Ring) Vector {
	a := MultiArray{
		coordMap: NewCoordMap([]int{len(elements)}),
		Array:    elements,
		baseRing: baseRing,
	}
	return Vector{&a}
}

// Populate initializes the items of this vector using a given function.
func (v Vector) Populate(f func(int) RingElement) Vector {
	for i := 0; i < v.Length(); i++ {
		v.SetElementAtIndex(i, f(i))
	}
	return v
}

// Map replaces every cell with the output of the given function in-place.
func (v Vector) Map(f func(RingElement, int) RingElement) Vector {
	return v.MultiArray.Map(func(el RingElement, coords []int) RingElement {
		return f(el, coords[0])
	}).AsVector()
}

// ForEach calls the given function with the contents of each cell.
func (v Vector) ForEach(f func(RingElement, int)) {
	v.MultiArray.ForEach(func(el RingElement, coords []int) {
		f(el, coords[0])
	})
}

// DotProduct performs a dot product of the vectors and returns the result.
func (v Vector) DotProduct(b Vector) Polynomial {
	out := NewPolynomial(v.baseRing)
	v.ForEach(func(el RingElement, i int) {
		el.MulAdd(b.ElementAtIndex(i), out)
	})
	return out
}
