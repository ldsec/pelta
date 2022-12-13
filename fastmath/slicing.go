package fastmath

// Slice represents a slice [start, end) over a list of values.
type Slice struct {
	start int
	end   int
}

func (p *Poly) Slice(s Slice) *Poly {
	return p

}

func (v *IntVec) Slice(s Slice) *IntVec {
	return v
}

func (v *PolyVec) Slice(s Slice) *PolyVec {
	return v
}
