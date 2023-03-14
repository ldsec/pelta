package logging

import (
	"fmt"
	"strings"
	"time"
)

var logging = true
var indentation = 0

type LogExecData struct {
	t0       time.Time
	procName string
	name     string
	short    bool
}

func getIndentation() string {
	return strings.Repeat("  ", indentation)
}

func indent() {
	indentation += 1
}

func unindent() {
	indentation -= 1
}

func IsProductionBuild() bool {
	return !logging
}

func Log(prefix, msg string) {
	if !logging {
		return
	}
	msg = strings.ReplaceAll(msg, "\n", fmt.Sprintf("\n%s%s: ", getIndentation(), prefix))
	fmt.Printf("%s%s: %s\n", getIndentation(), prefix, msg)
}

func LogWithoutNewline(prefix, msg string) {
	if !logging {
		return
	}
	msg = strings.ReplaceAll(msg, "\n", fmt.Sprintf("\n%s: ", prefix))
	fmt.Printf("%s%s: %s", getIndentation(), prefix, msg)
}

func PanicOnProduction(prefix, msg string) {
	Log(prefix, msg)
	if !logging {
		panic("execution failed!!")
	}
}

func LogExecStart(procName, description string) LogExecData {
	Log(procName, fmt.Sprintf("%s started", description))
	indent()
	return LogExecData{time.Now(), procName, description, false}
}

func LogExecShortStart(procName, description string) LogExecData {
	LogWithoutNewline(procName, fmt.Sprintf("%s...", description))
	indent()
	return LogExecData{time.Now(), procName, description, true}
}

func (l LogExecData) LogExecEnd() {
	if !logging {
		return
	}
	unindent()
	t1 := time.Now()
	if l.short {
		fmt.Printf("done (%dms)\n", t1.Sub(l.t0).Milliseconds())
		return
	}
	endStr := fmt.Sprintf("%s complete", l.name)
	dt := t1.Sub(l.t0).Milliseconds()
	updateData(l.procName, dt)
	Log(l.procName, fmt.Sprintf("%s (%dms)", endStr, dt))
}

func LogShortExecution(procName, description string, f func() interface{}) interface{} {
	if !logging {
		return f()
	}
	LogWithoutNewline(procName, fmt.Sprintf("%s...", description))
	t0 := time.Now()
	out := f()
	t1 := time.Now()
	dt := t1.Sub(t0).Milliseconds()
	updateData(procName, dt)
	fmt.Printf("done (%dms)\n", dt)
	return out
}

func LogExecution(procName, description string, f func() interface{}) interface{} {
	if !logging {
		return f()
	}
	Log(procName, fmt.Sprintf("%s started", description))
	indent()
	t0 := time.Now()
	out := f()
	t1 := time.Now()
	unindent()
	dt := t1.Sub(t0).Milliseconds()
	updateData(procName, dt)
	Log(procName, fmt.Sprintf("%s complete (%dms)", description, dt))
	return out
}
