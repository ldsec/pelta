package logging

import (
	"fmt"
	"strings"
	"time"
)

var logging bool = true
var indentation int = 0

type LogExecData struct {
	start  time.Time
	prefix string
	name   string
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
	return LogExecData{time.Now(), prefix, name}
}

func (l LogExecData) LogExecEnd() {
	unindent()
	Log(l.prefix, fmt.Sprintf("%s complete (%dms)", l.name, time.Now().Sub(l.start).Milliseconds()))
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
