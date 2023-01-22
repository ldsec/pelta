package logging

import (
	"fmt"
	"strings"
	"time"
)

var logging bool = true
var errorLogging bool = true
var indentation int = 0

type LogExecData struct {
	t0     time.Time
	prefix string
	name   string
	short  bool
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

func Log(prefix, msg string) {
	if !logging {
		return
	}
	msg = strings.ReplaceAll(msg, "\n", fmt.Sprintf("\n%s%s: ", getIndentation(), prefix))
	fmt.Printf("%s%s: %s\n", getIndentation(), prefix, msg)
}

// func LogError(prefix, msg string, err error, stopExecution bool) {
// 	if err != nil {
// 		if errorLogging {
// 			Log(fmt.Sprintf("!! %s", prefix), fmt.Sprintf("%s"))
// 		}
// 		if stopExecution {

// 		}
// 	}
// }

func LogWithoutNewline(prefix, msg string) {
	if !logging {
		return
	}
	msg = strings.ReplaceAll(msg, "\n", fmt.Sprintf("\n%s: ", prefix))
	fmt.Printf("%s%s: %s", getIndentation(), prefix, msg)
}

func LogExecStart(prefix, name string) LogExecData {
	Log(prefix, fmt.Sprintf("%s started", name))
	indent()
	return LogExecData{time.Now(), prefix, name, false}
}

func LogExecShortStart(prefix, name string) LogExecData {
	LogWithoutNewline(prefix, fmt.Sprintf("%s...", name))
	indent()
	return LogExecData{time.Now(), prefix, name, true}
}

func (l LogExecData) LogExecEnd() {
	unindent()
	t1 := time.Now()
	if l.short {
		fmt.Printf("done (%dms)\n", t1.Sub(l.t0).Milliseconds())
		return
	}
	endStr := fmt.Sprintf("%s complete", l.name)
	Log(l.prefix, fmt.Sprintf("%s (%dms)", endStr, t1.Sub(l.t0).Milliseconds()))
}

func LogShortExecution(prefix, name string, f func() interface{}) interface{} {
	if !logging {
		return f()
	}
	LogWithoutNewline(prefix, fmt.Sprintf("%s...", name))
	t0 := time.Now()
	out := f()
	t1 := time.Now()
	fmt.Printf("done (%dms)\n", t1.Sub(t0).Milliseconds())
	return out
}

func LogExecution(prefix, name string, f func() interface{}) interface{} {
	if !logging {
		return f()
	}
	Log(prefix, fmt.Sprintf("%s started", name))
	indent()
	t0 := time.Now()
	out := f()
	t1 := time.Now()
	unindent()
	Log(prefix, fmt.Sprintf("%s complete (%dms)", name, t1.Sub(t0).Milliseconds()))
	return out
}
