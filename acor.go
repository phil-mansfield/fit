package fit

/*
#cgo CFLAGS: -O2 -g
#cgo LDFLAGS: -lm

#include "acor.h"
*/
import "C"

import (
	"unsafe"
)

func AutocorrelationTime(x []float64) float64 {
	// TODO: implement in Go from first principles. Weare hasn't actually
	// released his code under any license, so we can't use it in a public
	// code.

	mean := 0.0
	sigma := 0.0
	tau := 0.0
	meanPtr := (*C.double)(unsafe.Pointer(&mean))
	sigmaPtr := (*C.double)(unsafe.Pointer(&sigma))
	tauPtr := (*C.double)(unsafe.Pointer(&tau))
	data := (*C.double)(unsafe.Pointer(&x[0]))

	C.acor(meanPtr, sigmaPtr, tauPtr, data, C.int(len(x)), C.int(10))

	return tau
}