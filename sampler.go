package fit

import (
	_ "fmt"
	"math"
	"math/rand"
	"time"
)

type Parameter struct {
	V float64 // Initial value.
	S float64 // Characteristic scale.
	Freeze bool
}

type SamplerConfig struct {
	Walkers int
	Stretch float64
	Rand *rand.Rand
}

var DefaultSamplerConfig = SamplerConfig{
	Walkers: 100,
	Stretch: 2.0,
	Rand: nil,
}

type Sampler struct {
	Walkers int // Number of walkers. Must be divisible by two.
	Chains  [][][][]float64
	Rand    *rand.Rand // Internal generator. Must be thread-safe.


	probs [][]float64 // ln(PDF) of the current location of each chain.
	nAccept, nSteps int
	// Buffered values.
	ap, api, asqrti, afact float64
}

func NewSampler(config ...SamplerConfig) *Sampler {
	instConfig := DefaultSamplerConfig
	if len(config) > 0 {
		instConfig = config[0]
	}

	rrand := rand.New(rand.NewSource(time.Now().UnixNano()))
	if instConfig.Rand != nil {
		rrand = instConfig.Rand
	}

	sampler := &Sampler{
		Walkers: instConfig.Walkers,
		Chains: [][][][]float64{
			make([][][]float64, instConfig.Walkers/2),
			make([][][]float64, instConfig.Walkers/2),
		},
		Rand: rrand,
		probs: [][]float64{
			make([]float64, instConfig.Walkers/2),
			make([]float64, instConfig.Walkers/2),
		},
		ap: instConfig.Stretch,
		api: 1.0/instConfig.Stretch,
		asqrti: 1.0/math.Sqrt(instConfig.Stretch),
		afact: instConfig.Stretch - 1,
	}

	return sampler
}

func (sampler *Sampler) Run(pdf LogPDF, param []Parameter, term ...Terminator) {
	if len(term) == 0 {
		panic("Default Terminator not implemented yet.")
	}

	sampler.init(pdf, param)

	for noStops(term, sampler.Chains) {
		sampler.step(pdf)
	}
}

func (sampler *Sampler) init(pdf LogPDF, param []Parameter) {
	p0 := make([]float64, len(param))
	for i := range p0 {
		p0[i] = param[i].V
	}

	for k := 0; k < 2; k++ {
		for j := range sampler.Chains[k] {
			pj := make([]float64, len(p0))
			copy(pj, p0)
			for i := range pj {
				pj[i] += (sampler.Rand.Float64()*2 - 1)*param[i].S
			}
			sampler.Chains[k][j] = append(sampler.Chains[k][j], pj)
			end := len(sampler.Chains[k][j]) - 1
			sampler.probs[k][j] = pdf(sampler.Chains[k][j][end])
		}
	}
}

func (sampler *Sampler) step(pdf LogPDF) {
	end := len(sampler.Chains[0][0]) - 1

	for k := 0; k < 2; k++ {
		kOther := 1 - k
		_ = kOther
		for j, chain := range sampler.Chains[k] {
			iTarget := sampler.Rand.Int63n(int64(len(sampler.Chains[0])))
			target := sampler.Chains[kOther][iTarget][end]
			lnp := sampler.probs[k][j]

			p, lnp := sampler.move(pdf, chain[end], target, lnp)

			sampler.probs[k][j] = lnp
			sampler.Chains[k][j] = append(chain, p)
		}
	}
}

func (sampler *Sampler) move(
	pdf LogPDF, src, target []float64, lnp float64,
) (p []float64, lnpNew float64) {

	dim := len(src)
	zf := sampler.afact * sampler.Rand.Float64()
	zr := (1 + zf)*(1 + zf)*sampler.api

	_ = zr
	p = make([]float64, dim)
	for i := range p {
		p[i] = target[i] + zr*(src[i] - target[i])
	}

	lnpNew = pdf(p)

	r := math.Log(sampler.Rand.Float64())
	if r < float64(dim-1)*math.Log(zr) + lnpNew - lnp {
		sampler.nSteps++
		sampler.nAccept++
		return p, lnpNew
	} else {
		copy(p, src)
		sampler.nSteps++
		return p, lnp
	}
}


func noStops(term []Terminator, chains [][][][]float64) bool {
	for _, t := range term {
		if t.Stop(chains) {
			return false
		}
	}
	return true
}