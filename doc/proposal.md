Proposal for a Go-Based Fitting Library
=======================================

In my opinion, one of the largest features missing from gonum is robust
curve-fitting support. This is a proposal designed to address that. It
will be broken up into three parts:

1. A general purpose background on fitting routines. This will likely be
review for most but will frame the argument that I make in the following
seciton.
2. An argument for structuring curve-fitting around Markov chain Monte Carlo
(with support from other methods where appropriate).
3. A proposed API for a curve-fitting package.

API Proposal
------------

The fundamental structure of this pacakge is as follows:
1. The core of the package is a `Sampler` struct, which implements the
Goodman-Weare MCMC sampling algorithm on a provided function pointer. This will
work very similarly to the battle-tested `EnsembleSampler` of the `emcee`
Python library, although Pythonic idoims will be replaced with Go idioms, and
several fundamental improvements (such as convergence testing and automatic
thinning/burning) will be added.
2. As a convenience feature, common PDFs will be provided. This will allow
users to fit arbitrary curves with 1D errors and upper/lower bounds, lines with
(potentially covaraint) 2D errors, multi-dimensional gaussian distributions,
curves with intrinsic scatter, and other common tasks without needing to 
write PDFs themselves. The only tasks which will not have default PDFs are
those which are sufficiently complicated that no existing fitting library
would have worked anyway.
3. Wrappers will be provided around all of these common PDFs so that novice
users can obtain curve fits with a single function call and will not need to
think about PDFs or distributions. In cases where a faster, case-specific
fitting routine is known, the wrapper can transparanetly switch the backend
to the faster method.

Package documentation will start by introducing the wrappers to keep from
overwhelming novice users.

Example Usages
--------------

Here is an example using the `Sampler` struct from my prototype of the
package's line-fitting wrapper. Note two things, first, I'm assuming a package
name of `fit` so this function call reads as `fit.Line(...)` in user code, and,
second, that this function is actually implementing four different fitting
routines, depending on what combination of x and y errors the user provides.

```go
func Line(x, y, xerr, yerr []float64) (param, err []float64, cov [][]float64) {
    // Pivoting code has been removed for clarity
    
    // PDF selection
    var pdf LogPDF
    switch {
    case xerr == nil && yerr == nil: pdf = LineConstantYErrorPDF(x, y)
    case xerr == nil: pdf = LineYErrorPDF(x, y, yerr)
    case yerr == nil: pdf = LineXErrorPDF(x, y yerr)
    default: pdf = LineXYErrorPDF(x, y, xerr, yerr)
    }
    
    // Calling the sampler
    sampler := NewSampler(pdf)
    y0, m0 := guessLineParameters(x, y)
    sampler.Run([]float64{y0, m0})
    return wrapperChainStats(sampler.Chains())
}
```

We'll get to the PDFs in a second, but let's unpack the last couple of lines of
the function. `sampler := NewSampler(pdf)` creates a `Sampler` for the given
PDF. The user may potentially want to provide configuration variables here,
but it's still an open quesiton what the best way to handle this is.
`y0, m0 := guessLineParameters(x, y)` obtains guesses for the line parameters.
This would usually need to be provided by the user, but we can get
good guesses without it for lines. `sampler.Run([]float64{y0, m0})` runs the
sampler until convergence with the given starting parameters. The user should
be able to specify exactly what convergence means (it will default to
Gelman-Rubin, but in principle the user could select something else like
some multiple of the autocorrelation time, a set number of steps, Raftery-Lewis,
or Geweke), but I'm not yet sure what the best way to specify this is.

So, to summarize, getting parameter distributions is just:
```go
sampler := NewSampler(pdf)
sampler.Run(paramGuess)
chains := sampler.Chains()
```

PDFs are represented by the `LogPDF` function type, `func([]float64) float64`.
These functions take a vector of parameters and return the logarithm
of the (not neccessarily normalized) probability of that parameter set. (The
logarithm is used because probabilities almost always involve evaluating an
exponential at the end, so this reduces the number of `math` calls needed and
preserves floating point accuracy.)

Open Questions in Package Design
--------------------------------

How to organize package-provided PDFs (structs vs. functions).

Line can take 2D errors, but Curve can't.

What should the defaults be?

Are we ever justified in switching to LM under the hood?

How to do configuration and convergence limit specification.

Which convergence estimators are worth implementing?
