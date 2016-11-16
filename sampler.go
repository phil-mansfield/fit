package fit

import (
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
	chains  [][]float64
	Rand    *rand.Rand // Internal generator. Must be thread-safe.

	skip, dim int
	probs  [][]float64 // ln(PDF) of the current location of each chain.
	nAccept, nSteps int
	// Buffered values.
	ap, api, asqrti, afact float64
	src, target, out []float64
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

// TODO: Stop people for passing chains that are too short.
func (sampler *Sampler) Run(pdf LogPDF, param []Parameter, term ...Terminator) {
	if len(term) == 0 {
		panic("Default Terminator not implemented yet.")
	}

	sampler.init(pdf, param)

	for noStops(term, sampler.chains) {
		sampler.step(pdf)
	}

	acor := sampler.AutocorrelationTimes()
	max := acor[0]
	for i := 1; i < len(acor); i++ {
		if acor[i] > max { max = acor[i] }
	}

	sampler.skip = int(math.Ceil(max))
	if sampler.skip < 0 { sampler.skip = 50 }
}

func (sampler *Sampler) init(pdf LogPDF, param []Parameter) {
	p0, pBuf := make([]float64, len(param)), make([]float64, len(param))
	for i := range p0 {
		p0[i] = param[i].V
	}

	sampler.dim = len(param)
	sampler.expandChains()

	for k := 0; k < 2; k++ {
		for j := 0; j < sampler.Walkers/2; j++ {
			copy(pBuf, p0)
			for i := range pBuf {
				pBuf[i] += (sampler.Rand.Float64()*2 - 1)*param[i].S
			}
			sampler.writeVec(k, 0, j, pBuf)
			sampler.probs[k][j] = pdf(pBuf)
		}
	}

	sampler.src = make([]float64, len(param))
	sampler.target = make([]float64, len(param))
	sampler.out = make([]float64, len(param))
}

func (sampler *Sampler) readVec(group, step, chain int, v []float64) {
	start := (group*sampler.Walkers/2 + chain)*sampler.dim
	copy(v, sampler.chains[step][start: start + sampler.dim])
}

func (sampler *Sampler) writeVec(group, step, chain int, v []float64) {
	start := (group*sampler.Walkers/2 + chain)*sampler.dim
	copy(sampler.chains[step][start: start + sampler.dim], v)
}

func (sampler *Sampler) step(pdf LogPDF) {
	end := len(sampler.chains) - 1
	sampler.expandChains()

	for k := 0; k < 2; k++ {
		kOther := 1 - k
		_ = kOther
		for j := 0; j < sampler.Walkers/2; j++ {
			iTarget := sampler.Rand.Int63n(int64(sampler.Walkers/2))
			sampler.readVec(kOther, end, int(iTarget), sampler.target)
			sampler.readVec(k, end, j, sampler.src)

			lnp := sampler.probs[k][j]

			lnpNew := sampler.move(
				pdf, sampler.src, sampler.target, sampler.out, lnp,
			)

			sampler.probs[k][j] = lnpNew
			sampler.writeVec(k, end+1, j, sampler.out)
		}
	}
}

func (sampler *Sampler) expandChains() {
	sampler.chains = append(
		sampler.chains, make([]float64, sampler.Walkers*sampler.dim),
	)
}

func (sampler *Sampler) move(
	pdf LogPDF, src, target, out []float64, lnp float64,
) (lnpNew float64) {

	dim := len(src)
	zf := sampler.afact * sampler.Rand.Float64()
	zr := (1 + zf)*(1 + zf)*sampler.api

	_ = zr
	for i := range out {
		out[i] = target[i] + zr*(src[i] - target[i])
	}

	lnpNew = pdf(out)

	r := math.Log(sampler.Rand.Float64())
	if r < float64(dim-1)*math.Log(zr) + lnpNew - lnp {
		sampler.nSteps++
		sampler.nAccept++
		return lnpNew
	} else {
		copy(out, src)
		sampler.nSteps++
		return lnp
	}
}


func noStops(term []Terminator, chains [][]float64) bool {
	for _, t := range term {
		if t.Stop(chains) {
			return false
		}
	}
	return true
}

func (sampler *Sampler) Chains() [][]float64 {
	nBurn := 20 * sampler.skip
	nIndependent := (len(sampler.chains)-nBurn) / sampler.skip + 1

	out := make([][]float64, sampler.dim)
	for i := range out {
		out[i] = make([]float64, 0, nIndependent * sampler.Walkers)
	}

	buf := make([]float64, sampler.dim)
	for group := 0; group < 1; group++ {
		for step := nBurn; step < len(sampler.chains); step += sampler.skip {
			for chain := 0; chain < sampler.Walkers/2; chain++ {
				sampler.readVec(group, step, chain, buf)
				for i := 0; i < sampler.dim; i++ {
					out[i] = append(out[i], buf[i])
				}
			}
		}
	}

	return out
}

// TODO: error checking
func (sampler *Sampler) AutocorrelationTimes() []float64 {


	avgChain := make([][]float64, sampler.dim)
	for i := range avgChain {
		avgChain[i] = make([]float64, len(sampler.chains))
	}
	buf := make([]float64, sampler.dim)

	for group := 0; group < 1; group++ {
		for step := 0; step < len(sampler.chains); step ++ {
			for chain := 0; chain < sampler.Walkers/2; chain++ {
				sampler.readVec(group, step, chain, buf)
				for i := 0; i < sampler.dim; i++ {
					avgChain[i][step] += buf[i]
				}
			}
		}
	}

	out := make([]float64, sampler.dim)
	for i := 0; i < sampler.dim; i++ {
		out[i] = AutocorrelationTime(avgChain[i])
	}

	return out
}