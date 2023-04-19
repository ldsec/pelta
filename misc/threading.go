package misc

import "sync"

var MaxRoutinesDefault = 16

// ExecuteInParallel executes the given `job` with the arguments 0, 1, ..., numCalls - 1
// in parallel and uses at most `maxRoutines` many goroutines
func ExecuteInParallel(numCalls int, job func(int), maxRoutines int) {
	doneChannel := make(chan struct{}, 2*maxRoutines)
	// initially, we can begin `maxRoutines` many jobs directly
	for i := 0; i < maxRoutines; i++ {
		doneChannel <- struct{}{}
	}
	var wg sync.WaitGroup
	wg.Add(numCalls)
	for i := 0; i < numCalls; i += 1 {
		// wait for a notification from the `doneChannel` before starting
		// a new job
		select {
		case <-doneChannel:
			go func(i int) {
				job(i)
				wg.Done()
				// indicate that a new job can be started
				doneChannel <- struct{}{}
			}(i)
		}
	}
	wg.Wait()
}
