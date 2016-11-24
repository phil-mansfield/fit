package fit

import (
	"math"
)

type LogPDF func(param []float64) float64
type Func1D func(param []float64, x float64) float64

func ConstantError1DPDF(x, y []float64, f Func1D) LogPDF {
	return func(param []float64) float64 {
		sum := 0.0
		p2, s := param[:len(param) - 1], param[len(param) - 1]
		if s < 0 { return math.Inf(-1) }
		for i := range x {
			dy := f(p2, x[i]) - y[i]
			sum += -dy*dy / (2*s*s) - math.Log(s)
		}
		return sum
	}
}

func Error1DPDF(x, y, yerr []float64, f Func1D) LogPDF {
	return func(param []float64) float64 {
		sum := 0.0
		for i := range x {
			dy := f(param, x[i]) - y[i]
			sy := yerr[i]
			sum += -dy*dy/(2*sy*sy)
		}
		return sum
	}
}

func ScatterError1DPDF(x, y, yerr []float64, f Func1D) LogPDF {
	return func(param []float64) float64 {
		sum := 0.0
		fp := param[:len(param) - 1]
		ms := param[len(param) - 1] // model sigma
		if ms < 0 { return math.Inf(-1) }

		for i := range x {
			dy := f(fp, x[i]) - y[i]
			ds := yerr[i] // data sigma
			//sNorm := ds*ds + ms*ms
			//sum += -0.5*(math.Log(sNorm) - dy*dy/sNorm)
			//s := math.Sqrt(ds*ds + ms*ms)
			s := ds*ds + ms*ms
			sum += -dy*dy / (2*s) - math.Log(math.Sqrt(s))
		}

		return sum
	}
}