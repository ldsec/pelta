package tests

import (
	"github.com/ldsec/codeBase/commitment/math"
	"math/big"
	"testing"
)

func TestPopulate(t *testing.T) {
	cm := math.NewCoordMap([]int{3, 5, 9, 7})
	m := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(1000))
		})
	for i := 0; i < len(m.Array); i++ {
		if m.Array[i].(*math.ModInt).Value.Cmp(big.NewInt(int64(i))) != 0 {
			t.Errorf("MultiArray.Populate")
		}
	}
}

func TestSetElements(t *testing.T) {
	cm := math.NewCoordMap([]int{6, 6})
	m := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(1000))
		})
	startCoords := []int{3, 4}
	endCoords := []int{4, 5}
	zero := math.NewModInt(0, big.NewInt(1))
	replacement := []math.RingElement{zero, zero, zero, zero, zero, zero, zero}
	m.SetElements(startCoords, endCoords, replacement)
	for i := 0; i < len(m.Array); i++ {
		currCoords := cm.ToCoords(i)
		x, y := currCoords[0], currCoords[1]
		if x >= 3 && x <= 4 && y >= 4 && y <= 5 && !(x == 4 && y == 5) {
			if m.Array[i].(*math.ModInt).Value.Cmp(&zero.Value) != 0 {
				t.Errorf("MultiArray.SetElements")
			}
		}
	}
}

func TestAdd(t *testing.T) {
	cm := math.NewCoordMap([]int{6, 6})
	m1 := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(2000))
		})
	m2 := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(2000))
		})
	m1.Add(m2)
	for i := 0; i < len(m1.Array); i++ {
		if m1.Array[i].(*math.ModInt).Value.Cmp(big.NewInt(int64(i*2))) != 0 {
			t.Errorf("MultiArray.Add")
		}
	}
}

func TestMul(t *testing.T) {
	cm := math.NewCoordMap([]int{6, 6})
	m1 := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(10000))
		})
	m1.Mul(math.NewModInt(10, big.NewInt(10000)))
	for i := 0; i < len(m1.Array); i++ {
		if !m1.Array[i].(*math.ModInt).Eq(math.NewModInt(int64(i*10), big.NewInt(10000))) {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestHadamard(t *testing.T) {
	cm := math.NewCoordMap([]int{6, 6})
	m1 := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(10000))
		})
	m2 := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(10000))
		})
	m1.Hadamard(m2)
	for i := 0; i < len(m1.Array); i++ {
		if m1.Array[i].(*math.ModInt).Value.Cmp(big.NewInt(int64(i*i))) != 0 {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestNeg(t *testing.T) {
	cm := math.NewCoordMap([]int{6, 6})
	m := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index), big.NewInt(1000))
		})
	m.Neg()
	for i := 0; i < len(m.Array); i++ {
		if m.Array[i].(*math.ModInt).Value.Cmp(big.NewInt(int64(-i))) != 0 {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestSum(t *testing.T) {
	cm := math.NewCoordMap([]int{2, 2})
	m := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index+1), big.NewInt(1000))
		})
	if m.Sum().(*math.ModInt).Value.Cmp(big.NewInt(10)) != 0 {
		t.Errorf("MultiArray.Sum")
	}
}

func TestProduct(t *testing.T) {
	cm := math.NewCoordMap([]int{2, 2})
	m := math.NewMultiArray(cm.Dims).
		Populate(func(dims []int) math.RingElement {
			index := cm.FromCoords(dims)
			return math.NewModInt(int64(index+1), big.NewInt(1000))
		})
	if m.Product().(*math.ModInt).Value.Cmp(big.NewInt(24)) != 0 {
		t.Errorf("MultiArray.Product")
	}
}

func TestForEach(t *testing.T) {
	cm := math.NewCoordMap([]int{2, 2})
	m := math.NewMultiArray(cm.Dims).Populate(func(dims []int) math.RingElement {
		index := cm.FromCoords(dims)
		return math.NewModInt(int64(index), big.NewInt(1000))
	})
	m.ForEach(func(el math.RingElement, coords []int) {
		i := cm.FromCoords(coords)
		if el.(*math.ModInt).Value.Cmp(big.NewInt(int64(i))) != 0 {
			t.Errorf("MultiArray.ForEach")
		}
	})
}

func TestMap(t *testing.T) {
	cm := math.NewCoordMap([]int{2, 2})
	m := math.NewMultiArray(cm.Dims).Populate(func(dims []int) math.RingElement {
		index := cm.FromCoords(dims)
		return math.NewModInt(int64(index), big.NewInt(1000))
	})
	m.Map(func(el math.RingElement, _ []int) math.RingElement {
		return el.Copy().Add(math.NewModInt(1, big.NewInt(1000)))
	})
	m.ForEach(func(el math.RingElement, coords []int) {
		i := cm.FromCoords(coords)
		if el.(*math.ModInt).Value.Cmp(big.NewInt(int64(i+1))) != 0 {
			t.Errorf("MultiArray.Map")
		}
	})
}
