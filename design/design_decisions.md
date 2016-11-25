### Organization of wrapper fitting functions

### Output format for fitting functions

(also, what is and isn't included and why)

### What should the output format report at the `Values` and `Error` parameters

### Which functional forms get their own wrapper functions/PDFs?

### p0 is a parameter to fit.Line.Foo()

### Fitting against X & Y errors for arbitrary curves is non-trivial

### Handling fits that insert parameters beyond what the user specifed

### Percentile functions

### Splitting `NewSampler` and `Sample` functions

### Optional parameters

(You can learn variadic options by example. And they're only something that
advanced users are using anyway.)

### Name of the `Samples` and `RawSamples`

### To what extent is it our responsibility to teach people how MCMC works?

(emcee doesn't do this, but there are certain tasks that require the user to
access the underlying sampler.)

### What should be the point where users are forced to 

(two conflicting things: it sucks if the user wants to do something simple and can
only do it by shifting to a more complicated API, but we also don't want to inflate our
public facing API by too much and we can't anticipate everything people might want to
do)

### How easy should it be to find the parameters of maximum likelihood?

### Built-in distribution fitting

### Is grouping fitting functions non-idiomatic?

### Should groups be named `FooGroup` or `Foo` or something else?

"funcs?"

### Should `WithScatter` be a model type or a method type?

### `Parameter` method signature

### How exactly should fititng power laws work?
