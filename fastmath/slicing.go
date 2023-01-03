package fastmath

// Slice represents a slice [start, end) over a list of values.
type Slice struct {
	start int
	end   int
}

func NewSlice(start, end int) Slice {
	return Slice{start, end}
}

func (s Slice) Size() int {
	return s.end - s.start
}

func (p *Poly) Slice(s Slice) *Poly {
	panic("not implemented")
}

func (v *IntVec) Slice(s Slice) *IntVec {
	sliced := NewIntVec(s.Size(), v.baseRing)
	for i := s.start; i < s.end; i++ {
		sliced.Set(i-s.start, v.Get(i))
	}
	return sliced
}

func (v *PolyVec) Slice(s Slice) *PolyVec {
	panic("not implemented")
}
