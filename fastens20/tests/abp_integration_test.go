package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestABPSimpleShortRing(t *testing.T) {
	bfvRing := fastmath.BFVFullShortCommtRing(7)
	m := bfvRing.D
	n := bfvRing.D
	tau := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithABP(tau, bfvRing.BaseRing.ModulusAtLevel[0], fastmath.NewSlice(0, n))
	params := fastens20.GeneratePublicParameters(config)
	if !fastens20.Execute(rel.S, params) {
		t.Errorf("execution failed!")
	}
}

func TestABPSimpleFullRing(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	tau := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithABP(tau, bfvRing.BaseRing.ModulusAtLevel[0], fastmath.NewSlice(0, n))
	params := fastens20.GeneratePublicParameters(config)
	if !fastens20.Execute(rel.S, params) {
		t.Errorf("execution failed!")
	}
}

func TestABPSimpleFullRingRebased(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	tau := 128
	rel := createRandomRelation(m, n, bfvRing)
	// Rebase
	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebased(commitmentRing)
	config := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel).
		WithABP(tau, commitmentRing.BaseRing.ModulusAtLevel[0], fastmath.NewSlice(0, n))
	params := fastens20.GeneratePublicParameters(config)
	if !fastens20.Execute(rebasedRel.S, params) {
		t.Errorf("execution failed!")
	}
}

func TestABPSlice(t *testing.T) {
	bfvRing := fastmath.BFVFullShortCommtRing(7)
	m := bfvRing.D
	n := bfvRing.D
	tau := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithABP(tau, bfvRing.BaseRing.ModulusAtLevel[0], fastmath.NewSlice(0, n/2))
	params := fastens20.GeneratePublicParameters(config)
	if !fastens20.Execute(rel.S, params) {
		t.Errorf("execution failed!")
	}
}
