package fastmath

// Slice represents a slice [start, end) over a list of values.
type Slice struct {
	Start int
	End   int
}

func NewSlice(start, end int) Slice {
	return Slice{start, end}
}

func (s Slice) Size() int {
	return s.End - s.Start
}

func (s Slice) Contains(index int) bool {
	return index >= s.Start && index < s.End
}

func (p *Poly) Slice(s Slice) *Poly {
	panic("not implemented")
}

func (v *IntVec) Slice(s Slice) *IntVec {
	// Optimization
	if s.Start%v.baseRing.N == 0 && s.End%v.baseRing.N == 0 {
		return v.SliceWithPolys(s.Start/v.baseRing.N, s.End/v.baseRing.N, s.Size())
	}
	sliced := NewIntVec(s.Size(), v.baseRing)
	for i := s.Start; i < s.End; i++ {
		sliced.SetCoeff(i-s.Start, v.GetCoeff(i))
	}
	return sliced
}

func (v *PolyVec) Slice(s Slice) *PolyVec {
	panic("not implemented")
}
