package math

import (
	"fmt"
	"strings"
)

type Vector struct {
	*MultiArray
}

// NewVectorFromSize constructs a new empty vector with the given size.
func NewVectorFromSize(dim int) Vector {
	return Vector{NewMultiArray([]int{dim})}
}

// NewVectorFromSlice constructs a new vector from the given slice.
// Warning: Does not copy the underlying elements.
func NewVectorFromSlice(elements []RingElement) Vector {
	a := MultiArray{
		coordMap: NewCoordMap([]int{len(elements)}),
		Array:    elements,
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

// Element returns the ith element of this vector.
func (v Vector) Element(i int) RingElement {
	return v.ElementAtIndex(i)
}

// Map replaces every cell with the output of the given function in-place.
func (v Vector) Map(f func(RingElement, int) RingElement) Vector {
	return v.MultiArray.Map(func(el RingElement, coords []int) RingElement {
		return f(el, coords[0])
	}).AsVec()
}

// ForEach calls the given function with the contents of each cell.
func (v Vector) ForEach(f func(RingElement, int)) {
	v.MultiArray.ForEach(func(el RingElement, coords []int) {
		f(el, coords[0])
	})
}

// Dot performs a dot product of the vectors and returns the result.
func (v Vector) Dot(b Vector) RingElement {
	if v.Length() != b.Length() {
		panic(fmt.Sprintf("Vector.Dot: Different sizes %d != %d", v.Length(), b.Length()))
	}
	out := v.ElementAtIndex(0).Copy().Zero()
	v.ForEach(func(el RingElement, i int) {
		el.MulAdd(b.ElementAtIndex(i), out)
	})
	return out
}

// Slice returns a view to the slice [start, end).
func (v Vector) Slice(start int, end int) Vector {
	if start < 0 || end > v.Length() || end < start {
		panic(fmt.Sprintf("Vector.Slice: Slice (%d, %d) out of bounds for %d", start, end, v.Length()))
	}
	return NewVectorFromSlice(v.Array[start:end])
}

// All returns true if all the elements return true on the given predicate.
func (v Vector) All(pred func(RingElement, int) bool) bool {
	return v.MultiArray.All(func(el RingElement, coords []int) bool {
		return pred(el, coords[0])
	})
}

// Any returns true if some element returns true on the given predicate.
func (v Vector) Any(pred func(RingElement, int) bool) bool {
	return v.MultiArray.Any(func(el RingElement, coords []int) bool {
		return pred(el, coords[0])
	})
}

// String returns a string representation of the vector.
func (v Vector) String() string {
	strs := make([]string, 0, v.Length())
	for _, e := range v.Array {
		strs = append(strs, e.String())
	}
	return fmt.Sprintf("(%d)[%s]", v.Length(), strings.Join(strs, ", "))
}
