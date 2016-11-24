package fit

import (
	"math"
)

func ConstantError1D(
	x, y []float64, p0 []Parameter, f Func1D,
) (p, err []float64, cov [][]float64)  {
	pdf := ConstantError1DPDF(x, y, f)
	sampler := NewSampler()
	sampler.Run(pdf, p0, Steps(20000))
	chains := sampler.Chains()
	return chainStats(chains)
}

func Error1D(
	x, y, yerr []float64, p0 []Parameter, f Func1D,
) (p, err []float64, cov [][]float64)  {
	pdf := Error1DPDF(x, y, yerr, f)
	sampler := NewSampler()
	sampler.Run(pdf, p0, Steps(20000))
	return chainStats(sampler.Chains())
}

func ScatterError1D(
	x, y, yerr []float64, p0 []Parameter, f Func1D,
) (p, err []float64, cov [][]float64)  {
	pdf := ScatterError1DPDF(x, y, yerr, f)
	sampler := NewSampler()
	sampler.Run(pdf, p0, Steps(20000))
	return chainStats(sampler.Chains())
}

func chainStats(chains [][]float64) (mean, err []float64, cov [][]float64) {
	dim := len(chains)

	mean = chainMean(chains)
	cov = chainCovariance(chains, mean)
	err = make([]float64, dim)
	for i := range err {
		err[i] = math.Sqrt(cov[i][i])
	}
	return mean, err, cov
}

func chainMean(chains [][]float64) []float64 {
	dim := len(chains)
	n := float64(len(chains[0]))

	means := make([]float64, dim)

	for c := range chains {
		offset := 0.0
		for i := range chains[c] {
			offset += chains[c][i]
		}
		offset /= n
		for i := range chains[c] {
			means[c] += chains[c][i] - offset
		}
		means[c] /= n
		means[c] += offset
	}

	return means
}

// TODO: rewrite to not run into floating point cancellation bugs in the wings
// of the covariance matrix.
func chainCovariance(chains [][]float64, mean []float64) [][]float64 {
	dim := len(chains)
	n := float64(len(chains[0]))

	cov := make([][]float64, dim)
	for i := range cov {
		cov[i] = make([]float64, dim)
	}

	for c1 := range chains {
		for c2 := c1 + 1; c2 < len(chains); c2++ {
			for i := range chains[0] {
				dx1 := chains[c1][i] - mean[c1]
				dx2 := chains[c2][i] - mean[c2]
				cov[c1][c2] += dx1*dx2
			}
			cov[c1][c2] /= n
		}
	}

	for c1 := range chains {
		for c2 := 0; c2 < c1; c2++ {
			cov[c1][c2] = cov[c2][c1]
		}
	}

	for c := range chains {
		for i := range chains[c] {
			dx := chains[c][i] - mean[c]
			cov[c][c] += dx*dx
		}
		cov[c][c] /= n
	}

	return cov
}
