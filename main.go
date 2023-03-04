package main

import (
	"fmt"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/ldsec/codeBase/commitment/relations"
)

func main() {
	logging.RegisterTimingCollection("Setup")
	logging.RegisterTimingCollection("ABPProver")
	logging.RegisterTimingCollection("ABPVerifier")
	for i := 0; i < 1; i++ {
		fmt.Println("**** RUN", i+1, "START")
		relations.RunKeyGenRelation()
		//relations.RunAggrKeyGenRelation()
		//relations.RunCollectiveDecRelation()
		//relations.RunAggrCollectiveDecRelation()
		//relations.RunKeySwitchRelation()
		//relations.RunCollectiveBootstrappingRelation()
		//relations.RunRelinKeyGenRelation()
		//relations.RunAggregationRelation()
		logging.FinishRun()
		fmt.Println("**** RUN", i+1, "END")
	}
	logging.FinishExperiment()
}
