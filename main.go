package main

import (
	"flag"
	"fmt"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/ldsec/codeBase/commitment/relations"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"strings"
)

func main() {
	// receive program arguments
	var numExecutions int
	var rlweDegree int
	var commtDegree int
	var levels int
	var delta1 uint64
	var lambda int
	var kappa int
	var relation string
	relationExecutors := map[string]func(){
		"keygen":         relations.RunKeyGenRelation,
		"aggr_keygen":    relations.RunAggrKeyGenRelation,
		"collective_dec": relations.RunCollectiveDecRelation,
		"keyswitch":      relations.RunKeySwitchRelation,
		"bootstrapping":  relations.RunCollectiveBootstrappingRelation,
		"relinkeygen":    relations.RunRelinKeyGenRelation,
	}
	var relationExecutorNames []string
	for k := range relationExecutors {
		relationExecutorNames = append(relationExecutorNames, k)
	}
	flag.IntVar(&numExecutions, "n", 1, "number of executions")
	flag.IntVar(&rlweDegree, "rlwe", 13, "log rlwe ring degree (10, 11, 12, 13, 14, or 15)")
	flag.IntVar(&commtDegree, "commt", 7, "commitment ring degree")
	flag.IntVar(&levels, "levels", 1, "number of levels for the ring (1, 2, or 3)")
	flag.Uint64Var(&delta1, "delta1", 25, "log delta1")
	flag.IntVar(&lambda, "lambda", 10, "lambda security parameter")
	flag.IntVar(&kappa, "kappa", 9, "kappa security parameter")
	flag.StringVar(&relation, "rel", "keygen", strings.Join(relationExecutorNames, ","))
	flag.Parse()
	logging.RegisterTimingCollection("Setup")
	logging.RegisterTimingCollection("ABPProver")
	logging.RegisterTimingCollection("ABPVerifier")
	// set the parameters
	var paramLiteral bfv.ParametersLiteral
	switch rlweDegree {
	case 10:
	case 11:
	case 12:
	case 13:
		paramLiteral = bfv.PN13QP218
		break
	case 14:
		paramLiteral = bfv.PN14QP438
		break
	case 15:
		paramLiteral = bfv.PN15QP880
		break
	default:
		panic("invalid rlwe ring degree, use --help to get the degrees")
	}
	if levels == 1 {
		relations.MainRing = fastmath.BFVZeroLevelRingCustom(paramLiteral, rlweDegree)
		relations.RebaseRing = fastmath.BFVZeroLevelRingCustom(paramLiteral, commtDegree)
	} else if levels == 2 {
		relations.MainRing = fastmath.BFVTwoLevelsRingCustom(paramLiteral, rlweDegree)
		relations.RebaseRing = fastmath.BFVTwoLevelsRingCustom(paramLiteral, commtDegree)
	} else if levels == 3 {
		relations.MainRing = fastmath.BFVFullRingCustom(paramLiteral, rlweDegree)
		relations.RebaseRing = fastmath.BFVFullRingCustom(paramLiteral, commtDegree)
	} else {
		panic("invalid ring level, use --help to get the available levels")
	}
	relations.Delta1 = uint64(1 << delta1)
	relations.Lambda = lambda
	relations.Kappa = kappa
	// get the executor
	executor, ok := relationExecutors[relation]
	if !ok {
		panic("invalid relation name, use --help to get the available relations")
	}
	for i := 0; i < numExecutions; i++ {
		fmt.Println("**** RUN", i+1, "START")
		executor()
		logging.FinishRun()
		fmt.Println("**** RUN", i+1, "END")
	}
	logging.FinishExperiment()
}
