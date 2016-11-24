# Proposal for a Go-Based Fitting Library

## Abstract

Curve-fitting is a common scientific task and a mainstay of numerical libraries. The current
absence of any Go-based general purpose curve fitting library is a major roadblock to the
acceptance of the language as a scientific programming platform.

This document will describe the structure of a proposed `fit` package to be added to gonum.
The core of the library will be a robust affine-invariant Markov chain Monte Carlo sampler,
although wrappers will be provided to hide this fact from non-expert users. The low-level
structure will be similar to the popular Python-based MCMC library
[emcee](http://dan.iel.fm/emcee/current/) (although lessons learned after a half decade of
user will be applied).

This document will be broken into two main parts. The first will list requirements and will
be accompanied by a code snippet that shows how this task can be accomplished in this library.
The second part will show the full API of the library.

## Requirements

* The most common fitting tasks must be trivial
* Placing limits on parameters
* Bounded but unknown points
* Freezing parameters
* Obtaining the full parameter distribution
* Finding the point of maximum likelihood
* Fitting a distribution
* Writing your own convergence checker
* Changing sampler parameters

### The most common fitting tasks must be trivial

Overwhelmingly, the most common fitting task is fitting an arbitrary curve to data
points with 1D errors.

```go
x, y, yerr := // ...

parabola := func(x, p []float64) float64 {
    return p[0] + x*p[1] + x*x*p[2]
}

// Initial guesses for parameters, can be very wrong.
p0 := []Parameter{ {V: 1}, {V: 2}, {V: 1.5}, }

out, err := fit.Curve.YErrors(parabola, p0, x, y, yerr)
if err != nil { panic(err.Error()) } // Could not converge
```

Here, `out` is a struct of type `fit.Output`

```go
type Output struct {
    Values, Errors []flaot64
    Covariance [][]float64
}
```

50% of people using this library won't need to know anything beyond this.

This can be generalized to other common tasks in the obvious way. For example:

```go
out, err = fit.Curve.UnknownError(parabola, p0, x, y)
out, err = fit.CurveWithScatter.XError(parabola, p0, x, y )
```

The latter case assumes that there is
some sort of intrinsic scatter to points in addition to known measurement errors
(something which is almost always true in, say, observed astrophysical relations).

There are some functional forms which can be fit more accurately if the fitting
routine knows what the function is at compile time.

```go
out, err = fit.Line.XYError(p0, x, y, xerr, yerr)
out, err = fit.PowerLaw.YError(p0, x, y, yerr)
```

### Placing limits on parameters

```go
x, y, yerr := // ...

parabola := func(x, p []float64) float64 {
    return p[0] + x*p[1] + x*x*p[2]
}

p0 := []Parameter{
    {V: 1}.UpperBound(3),
    {V: 2}.Limits(-5, +5),
    {V: 1.5}.Prior(myGaussian),
}

out, err := fit.Curve.YErrors(parabola, p0, x, y, yerr)
if err != nil { panic(err.Error()) } // Could not converge
```

### Freezing parameters

Fitting a curve to data is often an interactive process the requires a lot of
experimentation. It's often (read: almost always) useful to be able to "freeze"
parameters to a particular value and check how the fit behaves without allowing
that parameter to vary. Without library-level support for this, the user is
required to manually change the fitting function along with all analysis code,
which can be a massive pain.

```go
x, y, yerr := // ...

parabola := func(x, p []float64) float64 {
    return p[0] + x*p[1] + x*x*p[2]
}

p0 := []Parameter{
    {V: 1},
    {V: 2}.Freeze(), // Won't vary.
    {V: 1.5},
}

out, err := fit.Curve.YErrors(parabola, p0, x, y, yerr)
if err != nil { panic(err.Error()) } // Could not converge
```

### Obtaining the full parameter distribution

The `fit.Output` type has `Values`, `Errors`, and `Covariance` field. These are great
first order statistics of the fir, but don't tell the full story. What if you want to
know about 3-sigma error bars? What if parameter distribution isn't gaussian? What if you
want to use different statistics to represent errors? What if you want to check to see if
something went wrong? For this we will need to access the full poserior PDF. (At this point
we've far surpassed the functionality of most fitting libraries. If someone wants to go this
far, it's a safe bet that they know what MCMC is... and if they don't we should teach them.)

To do this we'll need to stop using the wrapper functions and directly call the underlying
sampler. This isn't too painful and consists of three steps: creating a sampler associated
without our data-curve combo, running it, and doing anlaysis of your choice on the resulting
samples:

```go
pdf := fit.CurvePDF.YErrors(parabola, x, y, yerr)
p0 := []Parameter{ {V: 1}, {V: 2}, {V: 1.5}, }

sampler := fit.NewSampler(pdf) // Make a new sampler
samples, err := sampler.Samples(p0) // Run until convergence
if err != nil { panic(err.Error()) } // Could not converge

// Do whatever analysis your heart desires on samples
```

Here we've replaced the `fit.Curve.YErrors` call with `fit.CurvePDF.YErrors` which
creates a probability distribution function that we can pass to `NewSampler`. The
resultant `samples` is a slice of parameter vectors which are each independently
drawn from this distribution. You can build histograms, take percentiles, etc.

### Finding the parameters of maximum likelihood

Nine times out of ten, you don't care about the most likely parameters, but about
properties of the global parameter distribution. However, sometimes you do (I won't
go into when here).

```go
pdf := fit.CurvePDF.YErrors(parabola, x, y, yerr)
p0 := []Parameter{ {V: 1}, {V: 2}, {V: 1.5}, }

sampler := fit.NewSampler(pdf) // Make a new sampler
_, err := sampler.Samples(p0) // Run until convergence
if err != nil { panic(err.Error()) } // Could not converge

p := sampler.MostLikelyParameters()
```

### Fitting a distribution



### Joint priors
### Writing your own convergence checker

`fit` will provide many built-in 

### Changing sampler parameters

```go
 sampler := NewSampler(pdf, fit.Walkers(150), fit.StepGranularity(1000))
```

## API

### Sampler

### Wrapper functions

```go
var (
	Curve Model
	Line
	PowerLaw
	CurveWithScatter
	LineWithScatter
	PowerLawWithScatter
)

type Model interface {
	UnknownErrors(x, y []float64) (Output, error)
	YErrors(x, y []float64) (Output, error)
	
	// Slow for arbitrary functions, but fast for power laws and lines:
	XYErrors(x, y []float64) (Output, error)
	// Can't be garuanteed to be correct for arbitrary curves:
	XErrors(x, y []float64) (Output, error)
}
```
