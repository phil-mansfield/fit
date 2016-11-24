package fit

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"
	"io/ioutil"
	"strings"
	"github.com/phil-mansfield/shellfish/math/sort"
	"github.com/phil-mansfield/shellfish/cmd/catalog"

	//plt "github.com/phil-mansfield/pyplot"
)

func NISTGaussianPDF(data []float64) LogPDF {
	return func(param []float64) float64 {
		mu, sigma := param[0], param[1]
		if sigma < 0 { return math.Inf(-1) }

		sum := 0.0
		for _, x := range data {
			dx := x - mu
			sum -= dx*dx
		}

		sum /= 2*sigma*sigma
		return sum - math.Log(sigma)*float64(len(data) + 1)
	}
}

///////////
// Tests //
///////////

func TestNIST(t *testing.T) {
	// This is the most difficult test in the NIST MCMC test suite.
	// See http://www.itl.nist.gov/div898/strd/mcmc/mcmc.html

	data := []float64{
		10000000000000.2,
		10000000000000.1,
		10000000000000.3,
		10000000000000.1,
		10000000000000.3,
		10000000000000.1,
		10000000000000.3,
		10000000000000.1,
		10000000000000.3,
		10000000000000.1,
		10000000000000.3,
	}

	pdf := NISTGaussianPDF(data)
	sampler := NewSampler()
	params := []Parameter{
		{V: 10000000000000, S: 1}, // mu
		{V: 1, S: 0.5}, // sigma
	}

	sampler.Run(pdf, params, Steps(20000))

	chains := sampler.Chains()
	mu, sigma := chains[0], chains[1]

	lead := 10000000000000.0
	mean := leadMean(mu, lead)
	std  := leadStd(mu, lead)
	p975 := sort.Percentile(mu, 0.975)
	p500 := sort.Percentile(mu, 0.5)
	p025 := sort.Percentile(mu, 0.025)

	// At 20000 steps, the accuracy is actually much higher, but I've
	// kept delta low so that every couple of tests doesn't get spurious
	// case failures (after all, there are ten values being checked).
	// Feel free to manually confirm this, but at 20000 steps we generate about
	// 2 million independent samples, so the error on the mean is about
	// one one thousandth of a standard deviation.
	delta := 0.005

	if !DeltaEq(mean, 10000000000000.2, delta) {
		t.Errorf("Got mean of mu %f, expected %f.",
			mean, 10000000000000.2)
	}
	if !DeltaEq(std, 0.033709993123162, delta) {
		t.Error("Got standard deviation of mu %f, expected %f.",
			std, 0.033709993123162)
	}
	if !DeltaEq(p975, 10000000000000.132819085883166, delta) {
		t.Errorf("Got 97.5%% of mu %f, expected %f.",
			p975, 10000000000000.132819085883166)
	}
	if !DeltaEq(p500, 10000000000000.2, delta) {
		t.Errorf("Got 50%% of mu %f, expected %f.",
			p500, 10000000000000.2)
	}
	if !DeltaEq(p025, 10000000000000.267180914116834, delta) {
		t.Errorf("Got 2.5%% of mu %f, expected %f.",
			p025, 10000000000000.26718091411683)
	}

	lead = 0.0
	mean = leadMean(sigma, lead)
	std  = leadStd(sigma, lead)
	p975 = sort.Percentile(sigma, 0.975)
	p500 = sort.Percentile(sigma, 0.5)
	p025 = sort.Percentile(sigma, 0.025)

	if !DeltaEq(mean, 0.108372230793914, delta) {
		t.Errorf("Got mean of mu %f, expected %f.",
			mean,  0.108372230793914)
	}
	if !DeltaEq(std, 0.027485625202104, delta) {
		t.Error("Got standard deviation of mu %f, expected %f.",
			std, 0.027485625202104)
	}
	if !DeltaEq(p975,  0.069871704416342, delta) {
		t.Errorf("Got 97.5%% of mu %f, expected %f.",
			p975, 0.069871704416342)
	}
	if !DeltaEq(p500, 0.103462818336964, delta) {
		t.Errorf("Got 50%% of mu %f, expected %f.",
			p500,  0.103462818336964)
	}
	if !DeltaEq(p025, 0.175493354741336, delta) {
		t.Errorf("Got 2.5%% of mu %f, expected %f.",
			p025, 0.1754933547413363)
	}
}

func DeltaEq(x, y, delta float64) bool {
	return x + delta > y && x - delta < y
}

func leadMean(xs []float64, lead float64) float64 {
	sum := 0.0
	for _, x := range xs { sum += (x - lead) }
	return sum / float64(len(xs)) + lead
}

func leadStd(xs []float64, lead float64) float64 {
	sum, sqrSum := 0.0, 0.0
	for _, x := range xs {
		sum += x - lead
		sqrSum += (x - lead)*(x -lead)
	}
	n := float64(len(xs))
	return math.Sqrt(sqrSum/n - sum*sum/(n*n))
}

func TestConstantError1D(t *testing.T) {
	f := func(p []float64, x float64) float64 {
		y0, m := p[0], p[1]
		return y0 + m*x
	}

	pTrue := []float64{3, -1}

	x := make([]float64, 30)
	for i := range x { x[i] = float64(i) / 15.0 - 1.0}
	y := make([]float64, 30)
	for i := range y { y[i] = f(pTrue, x[i])}

	err := 0.05
	rand.Seed(time.Now().UnixNano())

	sy := make([]float64, 30)
	for i := range y { sy[i] = y[i] + rand.NormFloat64()*err }

	p0 := []Parameter{{V: 2, S: 0.1}, {V: 0, S: 0.1}, {V:1, S:0.1}}
	p, std, cov := ConstantError1D(x, sy, p0, f)
	fmt.Println(p)
	fmt.Println(std)
	fmt.Println(cov)
}

func TestScatterError1D(t *testing.T) {
	f := func(p []float64, x float64) float64 {
		y0, m := p[0], p[1]
		return y0 + m*x
	}

	pTrue := []float64{3, -1, 0.3}

	x := make([]float64, 30)
	for i := range x { x[i] = float64(i) / 15.0 - 1.0}
	y := make([]float64, 30)
	for i := range y { y[i] = f(pTrue, x[i])}

	err := make([]float64, 30)
	for i := range err { err[i] = 0.1  + float64(i) / 100 }
	rand.Seed(time.Now().UnixNano())

	sy := make([]float64, 30)
	for i := range y {
		sy[i] = y[i] + rand.NormFloat64()*err[i] +
			rand.NormFloat64()*pTrue[2]
	}

	p0 := []Parameter{{V: 2, S: 0.1}, {V: 0, S: 0.1}, {V:1, S:0.1}}
	p, std, cov := ScatterError1D(x, sy, err, p0, f)
	fmt.Println(p)
	fmt.Println(std)
	fmt.Println(cov)
}

func TestNorrisNIST(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	txt, _ := ioutil.ReadFile("test_files/norris.dat")
	lines := strings.Split(string(txt), "\n")

	_, cols, _ := catalog.ParseCols(lines, []int{}, []int{0, 1})
	x, y := cols[0], cols[1]

	f := func(p []float64, x float64) float64 {
		return p[0] + p[1]*x
	}
	p0 := []Parameter{{V:0, S: 1}, {V:0, S:1}, {V:1, S:0.1}}

	p, std, cov := ConstantError1D(x, y, p0, f)
	fmt.Println(p)
	fmt.Println(std)
	fmt.Println(cov)
}

////////////////
// Benchmarks //
////////////////

func BenchmarkTenPointGaussian(b *testing.B) {
	pts := []float64{2, 1, 3, 1, 3, 3, 1, 3, 1, 3}
	pdf := NISTGaussianPDF(pts)
	sampler := NewSampler()
	params := []Parameter{ {V: 0.0, S: 1.0}, {V: 1.0, S: 0.5} }

	b.ResetTimer()
	steps := b.N / 100
	if steps < 1 { steps = 1 }
	sampler.Run(pdf, params, Steps(steps))
}

func BenchmarkTenPointGaussianPDF(b *testing.B) {
	pts := []float64{2, 1, 3, 1, 3, 3, 1, 3, 1, 3}
	pdf := NISTGaussianPDF(pts)
	p := []float64{1, 1}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		pdf(p)
	}
}