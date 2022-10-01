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
func NewVectorFromSlice(elements []Polynomial, baseRing *ring.Ring) Vector {
	a := MultiArray{
		coordMap: NewCoordMap([]int{len(elements)}),
		Array:    elements,
		baseRing: baseRing,
	}
	return Vector{&a}
}

// Populate initializes the items of this vector using a given function.
func (v Vector) Populate(f func(int) Polynomial) Vector {
	for i := 0; i < v.Length(); i++ {
		v.SetElementAtIndex(i, f(i))
	}
	return v
}

// Map replaces every cell with the output of the given function in-place.
func (v Vector) Map(f func(Polynomial, int) Polynomial) Vector {
	return v.MultiArray.Map(func(el Polynomial, coords []int) Polynomial {
		return f(el, coords[0])
	}).AsVector()
}

// ForEach calls the given function with the contents of each cell.
func (v Vector) ForEach(f func(Polynomial, int)) {
	v.MultiArray.ForEach(func(el Polynomial, coords []int) {
		f(el, coords[0])
	})
}
