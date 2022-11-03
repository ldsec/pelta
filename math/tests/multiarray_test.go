package tests

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"math/big"
	"testing"
)

func TestMultiArrayPopulate(t *testing.T) {
	cm := algebra.NewCoordMap([]int{3, 5, 9, 7})
	m := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(1000))
		})
	for i := 0; i < len(m.Array); i++ {
		if m.Array[i].(*rings.ZInt).Value.Cmp(big.NewInt(int64(i))) != 0 {
			t.Errorf("MultiArray.Populate")
		}
	}
}

func TestMultiArraySetElements(t *testing.T) {
	cm := algebra.NewCoordMap([]int{6, 6})
	m := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(1000))
		})
	startCoords := []int{3, 4}
	endCoords := []int{4, 5}
	zero := rings.NewZqInt(0, big.NewInt(1))
	replacement := []algebra.Element{zero, zero, zero, zero, zero, zero, zero}
	m.SetElements(startCoords, endCoords, replacement)
	for i := 0; i < len(m.Array); i++ {
		currCoords := cm.ToCoords(i)
		x, y := currCoords[0], currCoords[1]
		if x >= 3 && x <= 4 && y >= 4 && y <= 5 && !(x == 4 && y == 5) {
			if m.Array[i].(*rings.ZInt).Value.Cmp(zero.Value) != 0 {
				t.Errorf("MultiArray.SetElements")
			}
		}
	}
}

func TestMultiArrayAdd(t *testing.T) {
	cm := algebra.NewCoordMap([]int{6, 6})
	m1 := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(2000))
		})
	m2 := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(2000))
		})
	m1.Add(m2)
	for i := 0; i < len(m1.Array); i++ {
		if m1.Array[i].(*rings.ZInt).Value.Cmp(big.NewInt(int64(i*2))) != 0 {
			t.Errorf("MultiArray.Add")
		}
	}
}

func TestMultiArrayMul(t *testing.T) {
	cm := algebra.NewCoordMap([]int{6, 6})
	m1 := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(10000))
		})
	m1.Mul(rings.NewZqInt(10, big.NewInt(10000)))
	for i := 0; i < len(m1.Array); i++ {
		if !m1.Array[i].(*rings.ZInt).Eq(rings.NewZqInt(int64(i*10), big.NewInt(10000))) {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestMultiArrayHadamard(t *testing.T) {
	cm := algebra.NewCoordMap([]int{6, 6})
	m1 := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(10000))
		})
	m2 := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index), big.NewInt(10000))
		})
	m1.Hadamard(m2)
	for i := 0; i < len(m1.Array); i++ {
		if m1.Array[i].(*rings.ZInt).Value.Cmp(big.NewInt(int64(i*i))) != 0 {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestMultiArrayNeg(t *testing.T) {
	cm := algebra.NewCoordMap([]int{6, 6})
	m := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index+1), big.NewInt(1000))
		})
	m.Neg()
	for i := 0; i < len(m.Array); i++ {
		if m.Array[i].(*rings.ZInt).Value.Cmp(big.NewInt(int64(1000-i-1))) != 0 {
			t.Errorf("MultiArray.Hadamard")
		}
	}
}

func TestMultiArraySum(t *testing.T) {
	cm := algebra.NewCoordMap([]int{2, 2})
	m := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index+1), big.NewInt(1000))
		})
	if m.Sum().(*rings.ZInt).Value.Cmp(big.NewInt(10)) != 0 {
		t.Errorf("MultiArray.Sum")
	}
}

func TestMultiArrayProduct(t *testing.T) {
	cm := algebra.NewCoordMap([]int{2, 2})
	m := algebra.NewMultiArray(cm.Dims).
		Populate(func(dims []int) algebra.Element {
			index := cm.FromCoords(dims)
			return rings.NewZqInt(int64(index+1), big.NewInt(1000))
		})
	if m.Product().(*rings.ZInt).Value.Cmp(big.NewInt(24)) != 0 {
		t.Errorf("MultiArray.Product")
	}
}

func TestMultiArrayForEach(t *testing.T) {
	cm := algebra.NewCoordMap([]int{2, 2})
	m := algebra.NewMultiArray(cm.Dims).Populate(func(dims []int) algebra.Element {
		index := cm.FromCoords(dims)
		return rings.NewZqInt(int64(index), big.NewInt(1000))
	})
	m.ForEach(func(el algebra.Element, coords []int) {
		i := cm.FromCoords(coords)
		if el.(*rings.ZInt).Value.Cmp(big.NewInt(int64(i))) != 0 {
			t.Errorf("MultiArray.ForEach")
		}
	})
}

func TestMultiArrayMap(t *testing.T) {
	cm := algebra.NewCoordMap([]int{2, 2})
	m := algebra.NewMultiArray(cm.Dims).Populate(func(dims []int) algebra.Element {
		index := cm.FromCoords(dims)
		return rings.NewZqInt(int64(index), big.NewInt(1000))
	})
	m.Map(func(el algebra.Element, _ []int) algebra.Element {
		return el.Copy().Add(rings.NewZqInt(1, big.NewInt(1000)))
	})
	m.ForEach(func(el algebra.Element, coords []int) {
		i := cm.FromCoords(coords)
		if el.(*rings.ZInt).Value.Cmp(big.NewInt(int64(i+1))) != 0 {
			t.Errorf("MultiArray.Map")
		}
	})
}
