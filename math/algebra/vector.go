package algebra

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
func NewVectorFromSlice(elements []Element) Vector {
	a := MultiArray{
		coordMap: NewCoordMap([]int{len(elements)}),
		Array:    elements,
	}
	return Vector{&a}
}

// Populate initializes the items of this vector using a given function.
func (v Vector) Populate(f func(int) Element) Vector {
	for i := 0; i < v.Length(); i++ {
		v.SetElementAtIndex(i, f(i))
	}
	return v
}

// Element returns the ith element of this vector.
func (v Vector) Element(i int) Element {
	return v.ElementAtIndex(i)
}

// Map replaces every cell with the output of the given function in-place.
func (v Vector) Map(f func(Element, int) Element) Vector {
	return v.MultiArray.Map(func(el Element, coords []int) Element {
		return f(el, coords[0])
	}).AsVec()
}

// ForEach calls the given function with the contents of each cell.
func (v Vector) ForEach(f func(Element, int)) {
	v.MultiArray.ForEach(func(el Element, coords []int) {
		f(el, coords[0])
	})
}

// Dot performs a dot product of the vectors and returns the result.
func (v Vector) Dot(b Vector) Element {
	if v.Length() != b.Length() {
		panic(fmt.Sprintf("Vector.Dot: Different sizes %d != %d", v.Length(), b.Length()))
	}
	out := v.ElementAtIndex(0).Copy().Zero()
	v.ForEach(func(el Element, i int) {
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
func (v Vector) All(pred func(Element, int) bool) bool {
	return v.MultiArray.All(func(el Element, coords []int) bool {
		return pred(el, coords[0])
	})
}

// Any returns true if some element returns true on the given predicate.
func (v Vector) Any(pred func(Element, int) bool) bool {
	return v.MultiArray.Any(func(el Element, coords []int) bool {
		return pred(el, coords[0])
	})
}

// Diag converts this N-vector into a diagonal NxN matrix.
func (v Vector) Diag() Matrix {
	zero := v.Element(0).Copy().Zero()
	return NewMatrixFromDimensions(v.Length(), v.Length()).PopulateRows(
		func(i int) Vector {
			ithRow := NewVectorFromSize(v.Length()).Populate(
				func(_ int) Element {
					return zero.Copy()
				})
			ithRow.SetElementAtIndex(i, v.Element(i))
			return ithRow
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

func (v Vector) Appended(b Vector) Vector {
	appended := make([]Element, v.Length()+b.Length())
	for i := 0; i < v.Length(); i++ {
		appended[i] = v.MultiArray.Array[i]
	}
	for i := 0; i < b.Length(); i++ {
		appended[v.Length()+i] = b.MultiArray.Array[i]
	}
	return NewVectorFromSlice(appended)
}
