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

func (v *PolyVec) AddAll(q *Poly) *PolyVec {
	for _, p := range v.elems {
		p.Add(q)
	}
	return v
}

func (v *PolyVec) ScaleAll(factor uint64) *PolyVec {
	for _, p := range v.elems {
		p.Scale(factor)
	}
	return v
}

func (v *PolyVec) Add(b *PolyVec) *PolyVec {
	for i, p := range v.elems {
		p.Add(b.Get(i))
	}
	return v
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

func (v *PolyVec) Sum() *Poly {
	out := NewPoly(v.baseRing)
	for _, el := range v.elems {
		out.Add(el)
	}
	return out
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

func (v *PolyVec) NTT() *PolyNTTVec {
	nttPolys := make([]*PolyNTT, 0, v.Size())
	for _, p := range v.elems {
		nttPolys = append(nttPolys, p.NTT())
	}
	return &PolyNTTVec{nttPolys, v.baseRing}
}

// ToIntVec converts a polynomial vector of size m to an integer vector of size md where d is the polynomial degree.
func (v *PolyVec) ToIntVec() *IntVec {
	polys := v.elems
	return NewIntVecFromPolys(polys, v.Size()*v.baseRing.N, v.baseRing)
}

func (v *PolyVec) InfNorm() uint64 {
	max := v.elems[0].MaxCoeff(0)
	for _, p := range v.elems[1:] {
		c := p.MaxCoeff(0)
		if c > max {
			max = c
		}
	}
	return max
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

func (v *PolyNTTVec) AddAll(q *PolyNTT) *PolyNTTVec {
	for _, p := range v.elems {
		p.Add(q)
	}
	return v
}

func (v *PolyNTTVec) ScaleAll(factor uint64) *PolyNTTVec {
	for _, p := range v.elems {
		p.Scale(factor)
	}
	return v
}

func (v *PolyNTTVec) Add(b *PolyNTTVec) *PolyNTTVec {
	for i, p := range v.elems {
		p.Add(b.Get(i))
	}
	return v
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

func (v *PolyNTTVec) Sum() *PolyNTT {
	out := NewPoly(v.baseRing).NTT()
	for _, el := range v.elems {
		out.Add(el)
	}
	return out
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

func (v *PolyNTTVec) Eq(r *PolyNTTVec) bool {
	for i, p := range v.elems {
		if !p.Eq(r.Get(i)) {
			return false
		}
	}
	return true
}

func (v *PolyNTTVec) Mul(r *PolyNTTVec) *PolyNTTVec {
	for i, p := range v.elems {
		p.Mul(r.Get(i))
	}
	return v
}

func (v *PolyNTTVec) MulAll(r *PolyNTT) *PolyNTTVec {
	for _, p := range v.elems {
		p.Mul(r)
	}
	return v
}

func (v *PolyNTTVec) InvNTT() *PolyVec {
	polys := make([]*Poly, 0, v.Size())
	for _, p := range v.elems {
		polys = append(polys, p.InvNTT())
	}
	return &PolyVec{polys, v.baseRing}
}

func (v *PolyNTTVec) Dot(r *PolyNTTVec) *PolyNTT {
	out := NewPoly(v.baseRing).NTT()
	for i, p := range v.elems {
		a := p.actual.ref
		b := r.Get(i).actual.ref
		v.baseRing.MulCoeffsAndAdd(a, b, out.actual.ref)
	}
	return out
}

// ToIntVec converts a polynomial vector of size m to an integer vector of size md where d is the polynomial degree.
func (v *PolyNTTVec) ToIntVec() *IntVec {
	polys := make([]*Poly, 0, v.Size())
	for _, p := range v.elems {
		polys = append(polys, ForceInvNTT(p))
	}
	return NewIntVecFromPolys(polys, v.Size()*v.baseRing.N, v.baseRing)
}