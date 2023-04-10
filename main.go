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
	flag.IntVar(&rlweDegree, "rlwe", 13, "log rlwe ring degree")
	flag.IntVar(&commtDegree, "commt", 7, "commitment ring degree")
	flag.Uint64Var(&delta1, "delta1", 24, "log delta1")
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
	}
	relations.MainRing = fastmath.BFVZeroLevelRingCustom(paramLiteral, rlweDegree)
	relations.RebaseRing = fastmath.BFVZeroLevelRingCustom(paramLiteral, commtDegree)
	relations.Delta1 = delta1
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
