package tests

import (
	"github.com/ldsec/codeBase/commitment/math"
	"reflect"
	"testing"
)

func TestFromCoords(t *testing.T) {
	coordMap := math.NewCoordMap([]int{4, 5, 7})
	if coordMap.FromCoords([]int{0, 0, 0}) != 0 {
		t.Errorf("CoordMap.FromCoords: (0, 0, 0) => 0")
	}
	if !reflect.DeepEqual(139, coordMap.FromCoords([]int{3, 4, 6})) {
		t.Errorf("CoordMap.ToCoords: (3, 4, 6) => 139")
	}
	if coordMap.FromCoords([]int{3, 0, 0}) != 3 {
		t.Errorf("CoordMap.FromCoords: (3, 0, 0) => 3")
	}
	if coordMap.FromCoords([]int{3, 0, 0}) != 3 {
		t.Errorf("CoordMap.FromCoords: (3, 0, 0) => 3")
	}
	if coordMap.FromCoords([]int{0, 1, 0}) != 4 {
		t.Errorf("CoordMap.FromCoords: (0, 1, 0) => 4")
	}
	if coordMap.FromCoords([]int{1, 4, 0}) != 17 {
		t.Errorf("CoordMap.FromCoords: (1, 4, 0) => 17")
	}
	if coordMap.FromCoords([]int{1, 4, 3}) != 77 {
		t.Errorf("CoordMap.FromCoords: (1, 4, 3) => 77")
	}
}

func TestToCoords(t *testing.T) {
	coordMap := math.NewCoordMap([]int{4, 5, 7})
	if !reflect.DeepEqual([]int{0, 0, 0}, coordMap.ToCoords(0)) {
		t.Errorf("CoordMap.ToCoords: (0, 0, 0) <= 0")
	}
	if !reflect.DeepEqual([]int{3, 4, 6}, coordMap.ToCoords(139)) {
		t.Errorf("CoordMap.ToCoords: (3, 4, 6) <= 139")
	}
	if !reflect.DeepEqual([]int{3, 0, 0}, coordMap.ToCoords(3)) {
		t.Errorf("CoordMap.ToCoords: (3, 0, 0) <= 3")
	}
	if !reflect.DeepEqual([]int{3, 0, 0}, coordMap.ToCoords(3)) {
		t.Errorf("CoordMap.ToCoords: (3, 0, 0) <= 3")
	}
	if !reflect.DeepEqual([]int{0, 1, 0}, coordMap.ToCoords(4)) {
		t.Errorf("CoordMap.ToCoords: (0, 1, 0) <= 4")
	}
	if !reflect.DeepEqual([]int{1, 4, 0}, coordMap.ToCoords(17)) {
		t.Errorf("CoordMap.ToCoords: (1, 4, 0) <= 17")
	}
	if !reflect.DeepEqual([]int{1, 4, 3}, coordMap.ToCoords(77)) {
		t.Errorf("CoordMap.ToCoords: (1, 4, 3) <= 77")
	}
}

func TestReversed(t *testing.T) {
	coordMap := math.NewCoordMap([]int{4, 5, 7})
	rev := coordMap.Reversed()
	if !reflect.DeepEqual(rev.Dims, []int{7, 5, 4}) {
		t.Errorf("CoordMap.Reversed: (4, 5, 7) => (7, 5, 4)")
	}
	if !reflect.DeepEqual(rev.Mults, []int{1, 7, 35}) {
		t.Errorf("CoordMap.Reversed: (1, 4, 20) => (1, 7, 35)")
	}
}
