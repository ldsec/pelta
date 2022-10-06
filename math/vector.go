package math

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
func (v Vector) DotProduct(b Vector) RingElement {
	out := v.ElementAtIndex(0).Copy().Zero()
	v.ForEach(func(el RingElement, i int) {
		el.MulAdd(b.ElementAtIndex(i), out)
	})
	return out
}
