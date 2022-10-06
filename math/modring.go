package math

import (
	"math/big"
)

// ModInt represents an integer in the ring Z_q
// TODO utilize lattigo under the hood.
type ModInt struct {
	Value   big.Int
	Modulus big.Int
}

func NewModInt(value *big.Int, mod *big.Int) *ModInt {
	newValue, newMod := big.Int{}, big.Int{}
	newValue.Set(value)
	newMod.Set(mod)
	newValue.Mod(&newValue, &newMod)
	i := ModInt{newValue, newMod}
	return &i
}

func (m *ModInt) Add(q RingElement) RingElement {
	m.Value.Add(&m.Value, &q.(*ModInt).Value).Mod(&m.Value, &m.Modulus)
	return m
}

func (m *ModInt) Mul(q RingElement) RingElement {
	m.Value.Mul(&m.Value, &q.(*ModInt).Value).Mod(&m.Value, &m.Modulus)
	return m
}

func (m *ModInt) MulAdd(q RingElement, out RingElement) {
	// Keep m * q
	tmp := m.Copy()
	tmp.Mul(q)
	out.Add(tmp)
}

func (m *ModInt) Neg() RingElement {
	m.Value.Neg(&m.Value)
	return m
}

func (m *ModInt) Zero() RingElement {
	m.Value.Set(big.NewInt(0))
	return m
}

func (m *ModInt) One() RingElement {
	m.Value.Set(big.NewInt(1))
	return m
}

func (m *ModInt) Copy() RingElement {
	return NewModInt(&m.Value, &m.Modulus)
}
