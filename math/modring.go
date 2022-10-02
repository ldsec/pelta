package math

// ModInt represents an integer in the ring Z_q
// TODO utilize lattigo under the hood.
type ModInt struct {
	Value   int
	Modulus int
}

func NewModInt(value int, mod int) ModInt {
	return ModInt{value, mod}
}

func (m *ModInt) Add(q RingElement) RingElement {
	m.Value += q.(*ModInt).Value % m.Modulus
	m.Value = m.Value % m.Modulus
	return m
}

func (m *ModInt) Mul(q RingElement) RingElement {
	m.Value *= q.(*ModInt).Value
	m.Value = m.Value % m.Modulus
	return m
}

func (m *ModInt) MulAdd(q RingElement, out RingElement) RingElement {
	out.(*ModInt).Value += m.Value * q.(*ModInt).Value
	out.(*ModInt).Value = out.(*ModInt).Value % m.Modulus
	return m
}

func (m *ModInt) Neg() RingElement {
	m.Value = -m.Value % m.Modulus
	return m
}

func (m *ModInt) Zero() RingElement {
	m.Value = 0
	return m
}

func (m *ModInt) One() RingElement {
	m.Value = 1
	return m
}

func (m *ModInt) Copy() RingElement {
	new := ModInt{
		Value:   m.Value,
		Modulus: m.Modulus,
	}
	return &new
}
