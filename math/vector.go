package math

import (
	"github.com/ldsec/lattigo/v2/ring"
)

type Vector struct {
	*MultiArray
}

// NewVectorFromSize constructs a new empty vector with the given size.
func NewVectorFromSize(dim int, baseRing *ring.Ring) Vector {
	a := NewMultiArray([]int{dim}, baseRing)
	return Vector{&a}
}

// NewVectorFromSlice constructs a new vector from the given slice.
// Warning: Does not copy the underlying elements.
func NewVectorFromSlice(elements []Polynomial, baseRing *ring.Ring) Vector {
	a := MultiArray{
		coordMap: NewCoordMap([]int{len(elements)}),
		Array:    elements,
		baseRing: baseRing,
	}
	return Vector{&a}
}

// MapInPlace replaces every cell with the output of the given function in-place.
func (v Vector) MapInPlace(f func(Polynomial, int) Polynomial) {
	v.MultiArray.MapInPlace(func(el Polynomial, coords []int) Polynomial {
		return f(el, coords[0])
	})
}

// ForEach calls the given function with the contents of each cell.
func (v Vector) ForEach(f func(Polynomial, int)) {
	v.MultiArray.ForEach(func(el Polynomial, coords []int) {
		f(el, coords[0])
	})
}
