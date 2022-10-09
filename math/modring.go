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

func NewZeroModInt(mod *big.Int) *ModInt {
	newValue := big.NewInt(0)
	i := ModInt{*newValue, *mod}
	return &i
}

func NewOneModInt(mod *big.Int) *ModInt {
	newValue := big.NewInt(1)
	i := ModInt{*newValue, *mod}
	return &i
}

func NewModInt(value int64, mod *big.Int) *ModInt {
	i := ModInt{*big.NewInt(value), *mod}
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
	return NewModInt(m.Value.Int64(), &m.Modulus)
}

func (m *ModInt) Pow(exp int64) RingElement {
	out := m.Copy().One().(*ModInt)
	expAbs := exp
	if expAbs < 0 {
		expAbs *= -1
	}
	for i := int64(0); i < expAbs; i++ {
		out.Mul(m)
	}
	m.Value = out.Value
	// Try to invert if the exponent is negative.
	if exp < 0 {
		m.Inv()
	}
	return m
}

func (m *ModInt) Scale(factor uint64) RingElement {
	m.Value.Mul(&m.Value, big.NewInt(int64(factor)))
	return m
}

// Inv sets this integer to its multiplicative inverse and returns the result.
func (m *ModInt) Inv() *ModInt {
	m.Value.ModInverse(&m.Value, &m.Modulus)
	return m
}

// Uint64 returns the uint64 representation of this integer.
func (m *ModInt) Uint64() uint64 {
	return m.Value.Uint64()
}
