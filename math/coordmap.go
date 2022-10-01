package math

// CoordMap represents a mapping from a linear space to a multidimensional coordinate space.
type CoordMap struct {
	Dims  []int
	Mults []int
}

func NewCoordMap(dims []int) CoordMap {
	mults := []int{1}
	acc := 1
	// Calculate the base multiplicants.
	for i := 1; i < len(dims); i++ {
		acc *= dims[i-1]
		mults = append(mults, acc)
	}
	return CoordMap{dims, mults}
}

// -- Simple base conversion in-between the coordinate space and the index space.

func (m *CoordMap) FromCoords(coords []int) int {
	acc := 0
	for i := 0; i < len(m.Dims); i++ {
		acc += coords[i] * m.Mults[i]
	}
	return acc
}

func (m *CoordMap) ToCoords(linearIndex int) []int {
	acc := linearIndex
	coords := make([]int, len(m.Dims))
	for i := len(m.Mults) - 1; i >= 0; i-- {
		qtd := acc / m.Mults[i]
		rem := acc % m.Mults[i]
		coords[i] = qtd
		acc = rem
	}
	return coords
}

// Reversed returns a new coordinate map with reversed dimensions.
func (m *CoordMap) Reversed() CoordMap {
	// Reverse the dimensions.
	newDims := make([]int, len(m.Dims))
	for i := 0; i < len(newDims); i++ {
		newDims[i] = m.Dims[(len(newDims)-1)-i]
	}
	return NewCoordMap(newDims)
}

// Copied returns a copy of this coordinate map.
func (m *CoordMap) Copied() CoordMap {
	return CoordMap{
		Dims:  append([]int{}, m.Dims...),
		Mults: append([]int{}, m.Mults...),
	}
}
