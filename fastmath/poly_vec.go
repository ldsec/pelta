package fastmath

import "github.com/tuneinsight/lattigo/v4/ring"

type PolyVec struct {
	elems    []*Poly
	baseRing *ring.Ring
}

func NewPolyVec(size int, baseRing *ring.Ring) *PolyVec {
	elems := make([]*Poly, size)
	for i := 0; i < size; i++ {
		elems[i] = NewPoly(baseRing)
	}
	return &PolyVec{elems, baseRing}
}

func (v *PolyVec) Size() int {
	return len(v.elems)
}

func (v *PolyVec) Populate(f func(int) *Poly) {
	for i := range v.elems {
		v.elems[i] = f(i)
	}
}

func (v *PolyVec) All(pred func(int, *Poly) bool) bool {
	for i, p := range v.elems {
		if !pred(i, p) {
			return false
		}
	}
	return true
}

func (v *PolyVec) Get(i int) *Poly {
	return v.elems[i]
}

func (v *PolyVec) Set(i int, newElement *Poly) {
	v.elems[i] = newElement
}

func (v *PolyVec) Update(f func(int, *Poly) *Poly) {
	for i, vOld := range v.elems {
		v.elems[i] = f(i, vOld)
	}
}

func (v *PolyVec) Append(p *Poly) {
	v.elems = append(v.elems, p)
}

func (v *PolyVec) Copy() *PolyVec {
	newElems := make([]*Poly, 0, len(v.elems))
	for _, p := range v.elems {
		newElems = append(newElems, p.Copy())
	}
	return &PolyVec{newElems, v.baseRing}
}

func (v *PolyVec) Eq(r *PolyVec) bool {
	for i, p := range v.elems {
		if !p.Eq(r.Get(i)) {
			return false
		}
	}
	return true
}

// ToIntVec converts a polynomial vector of size m to an integer vector of size md where d is the polynomial degree.
func (v *PolyVec) ToIntVec() *IntVec {
	polys := v.elems
	return NewIntVecFromPolys(polys, v.Size()*v.baseRing.N, v.baseRing)
}

type PolyNTTVec struct {
	elems    []*PolyNTT
	baseRing *ring.Ring
}

func (v *PolyNTTVec) Size() int {
	return len(v.elems)
}

func (v *PolyNTTVec) Populate(f func(int) *PolyNTT) {
	for i := range v.elems {
		v.elems[i] = f(i)
	}
}

func (v *PolyNTTVec) All(pred func(int, *PolyNTT) bool) bool {
	for i, p := range v.elems {
		if !pred(i, p) {
			return false
		}
	}
	return true
}

func (v *PolyNTTVec) Get(i int) *PolyNTT {
	return v.elems[i]
}

func (v *PolyNTTVec) Set(i int, newElement *PolyNTT) {
	v.elems[i] = newElement
}

func (v *PolyNTTVec) Update(f func(int, *PolyNTT) *PolyNTT) {
	for i, vOld := range v.elems {
		v.elems[i] = f(i, vOld)
	}
}

func (v *PolyNTTVec) Append(p *PolyNTT) {
	v.elems = append(v.elems, p)
}

func (v *PolyNTTVec) Copy() *PolyNTTVec {
	newElems := make([]*PolyNTT, 0, len(v.elems))
	for _, p := range v.elems {
		newElems = append(newElems, p.Copy())
	}
	return &PolyNTTVec{newElems, v.baseRing}
}

// ToIntVec converts a polynomial vector of size m to an integer vector of size md where d is the polynomial degree.
func (v *PolyNTTVec) ToIntVec() *IntVec {
	polys := make([]*Poly, 0, v.Size())
	for _, p := range v.elems {
		polys = append(polys, ForceInvNTT(p))
	}
	return NewIntVecFromPolys(polys, v.Size()*v.baseRing.N, v.baseRing)
}