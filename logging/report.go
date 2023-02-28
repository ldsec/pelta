package logging

import (
	"fmt"
	"math"
	"strings"
)

type ProcTiming struct {
	procName string
	time     int64
}

var registeredPrefixes = make(map[string][]ProcTiming)

func RegisterTimingCollection(procNamePrefix string) {
	registeredPrefixes[procNamePrefix] = []ProcTiming{}
}

func updateData(procName string, dt int64) {
	for k, v := range registeredPrefixes {
		if strings.HasPrefix(procName, k) {
			registeredPrefixes[k] = append(v, ProcTiming{procName: procName, time: dt})
		}
	}
}

type ExperimentStats struct {
	procNamePrefix string
	procStats      []ProcStatistics
	avg            float64
	std            float64
}

func (s ExperimentStats) String() string {
	outStrBuilder := strings.Builder{}
	outStrBuilder.WriteString(fmt.Sprintf("Report for %s\n\t", s.procNamePrefix))
	for i, procStat := range s.procStats {
		outStrBuilder.WriteString(fmt.Sprintf("Run %d: %v\n\t", i+1, procStat))
	}
	outStrBuilder.WriteString(fmt.Sprintf("Avg = %v\n\t", s.avg))
	outStrBuilder.WriteString(fmt.Sprintf("Std = %v\n\t", s.std))
	return outStrBuilder.String()
}

type RunStats struct {
	stats []ProcStatistics
}

type ProcStatistics struct {
	procNamePrefix string
	totalTime      int64
	timings        []ProcTiming
}

var runStats []RunStats

func computeProcStatistics(procNamePrefix string, timings []ProcTiming) ProcStatistics {
	sum := int64(0)
	for _, t := range timings {
		sum += t.time
	}
	return ProcStatistics{procNamePrefix: procNamePrefix, totalTime: sum, timings: timings}
}

func computeExperimentStatistics(procNamePrefix string, statsAcrossRuns []ProcStatistics) ExperimentStats {
	sum := int64(0)
	for _, s := range statsAcrossRuns {
		sum += s.totalTime
	}
	avg := float64(sum) / float64(len(statsAcrossRuns))
	tmp := float64(0)
	for _, s := range statsAcrossRuns {
		tmp += math.Pow(float64(s.totalTime)-avg, 2)
	}
	std := math.Sqrt(tmp / float64(len(statsAcrossRuns)))
	return ExperimentStats{procNamePrefix: procNamePrefix, procStats: statsAcrossRuns, avg: avg, std: std}
}

func FinishRun() {
	// collect the process stats
	var stats []ProcStatistics
	for prefix, timings := range registeredPrefixes {
		stats = append(stats, computeProcStatistics(prefix, timings))
	}
	// clear the timings
	for prefix := range registeredPrefixes {
		registeredPrefixes[prefix] = nil
	}
	// save the run stats
	runStats = append(runStats, RunStats{stats: stats})
}

func FinishExperiment() {
	for k := range registeredPrefixes {
		var procStatsAcrossRuns []ProcStatistics
		for _, run := range runStats {
			for _, procStat := range run.stats {
				if procStat.procNamePrefix == k {
					procStatsAcrossRuns = append(procStatsAcrossRuns, procStat)
				}
			}
		}
		// compute results across runs
		expStats := computeExperimentStatistics(k, procStatsAcrossRuns)
		fmt.Println(expStats)
	}
}
