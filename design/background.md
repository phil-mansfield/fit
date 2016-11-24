# Background

This section will have two subsections: the first focusing on the most common fitting
approach found in numerical libraries, least squares, and a second focusing
on the proposed approach for the `fit` pacakge, Markov chain Monte Carlo.

This can be safely skipped if you are already familiar with the motivation behind using
Markov chains to fit curves or if you're willing to trust me.

## Least-Squares

The most well-known fitting approach is least-squares minimization. Here, the best-fit
parameters of a function are taken to be the ones which result in the smallest squared
residuals. Formally, suppose we're fitting the 2D function `f(x; p)`, where `p` are the
unknown parameters: least squares methods would attempt to pick out the value of
`p` that minimizes the function `sum_i (f(x_i; p) - y_i)^2` across our data set
`{(x_i, y_i)`. The simplest least squares algorithm would be to create big gird of `p`
vectors over the entire viable parameter range, to evaluate the squared residual sum
for every vector and return the vector that gives the smallest sum. There are a number of
[much more sophisticated searching algorithms](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm)
than this. Most will merely find a local minimum (as opposed to a global minimum), but are
very fast.

The choice to minimize the squared residuals (as opposed to, say, residuals to the fourth
power) is not arbitrary: it can be understood in rigorous Bayesian terms as the assumption
of a very specific (and fairly reasonable) model for the errors associated with each point.
It assumes perfect accuracy in the `x` coordinates of every point and unknown but
constant-width gaussian errors on `y` coordiantes of every point (i.e. that
`P(y_true | y_measured)` is a guassian centered on `y_measured`), and no prior information
on what the reasonable value range for our parameters is. (There is an alternate version
of least squares called weighted least squares that can be used if there are known
inhomogenous gaussian error bars on every `y` measurement. The derivation is the same) Bayes'
theorem tells us that the probability of a parameter vector given the data we measured is
given by `P(p | data) = P(data | p) * P(p) / P(data)`. Since we assumed no priors on
`p`, the maximum of `P(p | data)` will be the maximum of `P(data | p)`. Expanding ou
that term gives us that
`P(data | p) = prod_i P((x_i, y_i) | p) = prod_i {exp(-(f(x; p) - x_i)^2/(2*sigma^2)) / (sqrt(2) * sigma)} `
due to our assumption about measurement errors (here `sigma` is that unknown but constant
error I mentioned earlier). The location of the maximum of `P(data | p)` is the same as
the location of the maximum of `ln(P(data | p))` and all constant offsets and constant
scalings have no effect on maximum's location, so the maximum of `P(p | data)` is the
same as the maximum of `-sum_i (f(x_i; p) - y_i)^2`. Ergo, the most likely parameter set
is the one that minimizes the squared residuals. (You may have noticed a small lie in this
derivation. I'll get to that later.)

There are two classes of issues associated with least squares minimization. The first is
issues that arise due to applying least squares to data sets where the underlying
statistical mode is not correct (I'll call these issues "class A") and the second is
problems that arise even when this model is correct ("class B"). I'll list the most
important class A and class B issues.

###Class A:

* _Non-gaussian errors_ - The central limit theorem means that the assumption of guassian
error bars are usually a pretty good one, but there are two common cases where this
isn't true. The first is when you've log-scaled your data before fitting (a common
procedure if you're fitting power laws or exponentials): you can use error propagation to
rescale the error bars on your points, but the assumption of a gaussian funcitonal form
for `P(x_i | p)` is no longer true. The second is when points have asymmetric error bars.
* _2D Error bars_ - There isn't a great way to fit curves to data with 2D (or, for that
matter, ND)
errors with least squares. There's a potpourri of ways to attempt doing this, but I'm not
going to discuss them in detail here: in general their derivations (if they aren't purely heuristic) make
very rescrictive assumptions about the relative size of x and y errors and even then are
only returning approximations to the most likely parameter vector.
* _Covariant error bars_ - covariant errors are even more difficult to handle than 2D
errors. (I don't actually know of any popular least squares fitting libraries that
allow the user to specify covariances.)
* _Priors_ - Priors can be thought of as a more general version of parameter bounds. Simple
1D cutoffs can be implemented can be hacked into least squares methods without any issues,
but more complicated cutoff shapes are more difficult to deal with. If you have priors that
aren't just sharp cutoffs (for example, if you have other measurements of a parameter), there's
no way to handle them.
* _Upper and lower limits_ - There is no robust, fully consistent way to handle limts with
least squares. They essentially have to be treated as very complicated priors on the
parameter space, which least squares is already not very good at dealing with. If your
upper limits aren't sharp cutoffs, it becomes impossible to model them.
* _Distributions_ - The entire set up described in the derivation above assumes that the
true location of each point is somewhere along an infinitely thin line and is then
scattered from its true location according to its measured error. However, _many_ physical
relations have intrinsic scatter (such as just about any observation relation in
astronomy/astrophysics). If your statistical model doesn't take this into account (and least
squares cannot), you won't be maximizing the correct probability function. In practice this
can lead to overly wide parameter distributions or artificial local maxima.
With sufficient [elbow grease](http://scipy-cookbook.readthedocs.io/items/FittingData.html#fitting-a-2d-gaussian),
least squares can be coaxed into fitting very simple distributions to data, but it is not
an easy or natural task.

###Class B:

* _Local Minima are not Global Minima_ - This is an issue which is brought up frequently (maybe
too frequently), but all popular least squares finders stop upon finding a local minimum, which
might not be a global minimum. Additionally, there will be nothing about the output of the finder
which could alert the user that this has happened.
* _Numerical Derivatives are Messy_ - The best least squares estimators rely on knowing the
Jacobian of the fitted function. If this can't be provided analytically it will need to be
estimated numerically. This is an issue because the value the objective function takes on in the
outskirts can be noisy, which can sometimes cause the estimator to aimlessly wander around the outskirts of
parameter space wihtout converging.
* _Unknown Error Bars are a Model Parameter_ - You may have noticed that in the earlier derivation
I was playing fast and loose with the error term. In reality what I should have done is computed
`P(p, sigma | data)` instead of `P(p | data)`. The reason I didn't do this is because if you had,
I would have found that the point of maximum likelihood didn't neccessaily need to coexist with
the point that minimizes square residuals.
* _Points of Maximum Likelihood are Less Useful Than You Might Think_ - While it might sound good
to pick out the point of maximum likelihood, this isn't always the way you want to report
best-fit parameters. There are a few toubling aspects to using points of maximum likelihood
(e.g. What if the marginalized point of maximum likelihood for a given parameter is in a
different location than the multivariate point of maximum likelihood? What if the objective
function is bimodal and the integrated likelihood under the lower peak is higher than the
integrated likelihood under the higher peak?) and there are other ways of reporting
parameters which have attractive propreties (such as taking the median of the marginalized
PDF for each parameter). With least squares you do not have the option to 
* _Least Squares Methods Aren't Good at Reporting Errors_ - Least squares methods focus on finding
the peak of `P(p | data)`, but errors on fit parameters are determined by the wings of the
`P(p | data)`, where most least squares estimators don't have much information. The standard approach
is to evaluate the Jacobian near the peak, to model `P(p | data)` as a perfect Gaussian and to
extrapolate out to the 68% contours. This works reasonably well if you're only interested in
1-sigma errors, but if you wanted to make statements like "This fit rules out theory X to an accuracy
of 3 sigma," you won't be able to. Furthermore, you have no indication of how accurate your
error estimates.

(These are not obscure problems. Every "Class B" issue I described here show up in many of the NIST
least squares [correctness tests](http://www.itl.nist.gov/div898/strd/general/dataarchive.html).)

None of this is to say that least squares is "wrong." Least squares estimators do exactly
they are supposed to do (finding the point of maximum likelihood for the simplest error model),
and they do it very quickly. It's just that they should not be viewed
as "general purpose" fitters. While our `fit` package should certainly contain some sort of
Levenberg-Marquardt functionality, the core behavior will need to be provided by something else.

## Markov Chains

An obvious alternative to finding just the point of maximum likelihood is to evaluate the entire
PDF. Naively, one could imagine making a grid across the viable parameter space
