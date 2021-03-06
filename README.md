# WeibullCount.jl

`WeibullCount.jl` provides a (partial) implementation of `WeibullCountModel` as a `Distribution` (in `Distributions.jl`).

A **Weibull count model** is a discrete distribution obtained as the *count model with inter-arrival times assumed to be i.i.d. Weibull distributed* [1, Section 3.1]. In words, a Weibull count model models the *number of events* which occur within some time `t`, where the time *between* each of the events is (i.i.d.) Weibull distributed.

## Examples
```julia
# these specific parameters caused issues with my Python implentation
x = 2
c = 0.36787944117144233
λ = 20.08553692

d = WeibullCountModel(λ, c) 
# => WeibullCountModel(λ=2.0085536919999999128094714251346886157989501953125e+01, c=3.67879441171442334024277442949824035167694091796875e-01, t=1.0)

pdf(d, x)
# => 1.953731363534023055349256598763820641174059397326485646480541542361423252806044e-09

pdf(d, [0, 1, 2])
# => 3-element Array{BigFloat,1}:
#   1.521522983183341907697464611074545778903924912431785205683339636283913506010731e-09
#   1.844778178908130394109566390808204884626471864468498910516595586139864555257461e-10
#   1.953731363534023055349256598763820641174059397326485646480541542361423252806044e-09

# the following required upping `N_TERMS` from 300 to 500, otherwise it did not converge
@time pdf(d, 20)
# => 138.153275 seconds (261.55 M allocations: 8.622 GiB, 1.67% gc time)
#    5.116056963580013944023420716711374369391236996632361437953045961338664334656706e-04

# after computing the case of `x = 20`
@time pdf(d, 15)
# => 0.225187 seconds (4.43 M allocations: 270.317 MiB, 14.76% gc time)
#    2.519010611989582537366569841564264788689589089612177323202223626875218163101008e-09
```

## Numerical stability
All computations are performed using `BigInt` and `BigFloat`, and will throw a `DivergedError` upon convergence-failure. The default number of terms in used for the approximation is given by `WeibullCount.N_TERMS = 300`, which I've found to work fairly well. 

The series is approximated by incrementally including more terms, and returning early if the result is a valid probability and the difference between the current and previous iteration is sufficiently small. see `WeibullCount.weibull_count_pdf_approx` for more info and implementation.

## Efficiency
In general, computing `pdf` of a Weibull count model is *slow*. I've tried to take measures to make it as efficient, but still, one can only do so much. It's mostly useful for small values of `x`.

### Euler (van Wijngaarden) transform
We make use of series acceleration by Euler transforming the alternating series. In particular, we use the van Wijngaarden transform, which is an efficient numerical implementation of the Euler transform.

### Caching
Results from `WeibullCount.weibull_alpha` are cached using an `LRUCache` from `LRUCache.jl`, which has default size of `10000`. This cache can be cleared and resized using `WeibullCount.cache_clear!()` and `WeibullCount.cache_size!(n)`, respectively.

This computes the `αⱼⁱ` in Eq. (11) in [1]. Due to its recursive nature, and having two of three arguments being integers, caching the results can speed up the computation immensely. As a result, if one is interested in optimizing an objective involving `pdf` of `WeibullCountModel` wrt. parameters `λ` and `c`, e.g. MLE estimate, it is advisable to formulate the computation such that successive calls to `pdf` (and thus `WeibullCount.weibull_alpha`) are performed with *fixed* `c`. An example might be to use take a round-robin approach, first optimizing wrt. `λ` while keeping `c` fixed, and then optimize wrt. `c` while keeping `λ` fixed.

## Motivation
The PDF of a Weibull count model is a given by an alternating series, and can, depending on choice of parameters `λ, c, t`, be very difficult to approximate numerically. The terms involve multiple (and increasing) factors of the `Γ`, i.e the Gamma function, and can therefore have very large absolute values.

First I obtained an implementation in Python, making extensive use of `log` for numerical stability and series acceleration (specifically, the van Wijngaarden transform, an implementation of Euler transform) to reduce number of necessary terms required in the partial sum to obtain covergence. This was to no avail. And this was not just a matter of increasing the number of terms in the partial sum, as this only resulted in `NaN`. This all pointed towards numerical instabilities being the main source of the issue. Furthermore, I was interested in performing an MLE estimate of the parameters of a Weibull count model using `scipy.optimize.minimize`, which in every case lead to the optimizer simply "exploiting" the numerical inaccuracies and making the partial sum blow up. 

Unfortunately, at least not that as far as I'm aware, one cannot simply replace all ones computations in Python with higher-precision integers and floats. Some methods do indeed provide a good flexibility, but most do not. And, as I found out later with this implementation, even 128-bit floats where insufficient for convergence in some cases, and so `np.float128` would not even be cut out for the job.

I therefore reached for Julia with its fantastic support for arbitrary precision through `BigInt` and `BigFloat` and parametric types; and here we are.

## Some maths

The Weibull counting model is a probability distribution derived as a /count model/ where the interarrival times are assumed to be i.i.d. Weibull distributed, i.e. interarrival times have PDF

    f(t) = λ c t^(c - 1) e^(- λ t^c)

and CDF

    F(t) = 1 - e^(- λ t^c)

with c, λ ∈ ℝ⁺.

A *count model* is obtained as follows. Let
- `Tᵢ` is the time at which the i-th even occurs
- `X(t)` is the number of events that occurred *up until time `t`*.

Then

    Tᵢ ≤ t    ⟺    X(t) ≥ i

i.e. `X(t) ≥ i` can only happen if the time at which the i-th event occurred is ≤ t, and vice versa.

This gives us the relationship

    Cᵢ(t) = P{ X(t) = i } = P { X(t) ≥ i } - P { X(t) ≥ i + 1 } = P { X(t) ≤ i  } - P { X(t) ≤ i + 1 }

where `Cᵢ(t)` then denotes the probability of `i` events having occurred at time `t`. Letting `Fᵢ(t)` denote the CDF of `Tᵢ`, we can write

    Cᵢ(t) = P{ X(t) = i } = F_{i}(t) - F_{i + 1}(t)
    
Therefore we need an expression for `F_{i}(t)`. Observe that we have

    F₁(t) = ∫ f(s) d(s),                     s ∈ [0, t]
    F₂(t) = ∫ ∫ f(s₁) f(s₂) d(s₁) d(s₂),     s₂ ∈ [0, t - s₁],    s₁ ∈ [0, t]

where `f(s)` denote the PDF of a Weibull random variable. 
But note that the inner-most intergral in `F₂(t)` is just `F₁(t - s₁)`, and so we have

    F₂(t) = ∫ F₁(t - s) f(s) d(s),    s ∈ [0, t]
    
For any general `i` we have

    Fᵢ(t) = ∫ F_{i - 1}(t - s) f(s) d(s),    s ∈ [0, t]

In [1] they substitute this into the expression for `Cᵢ(t)` above, and consider the Taylor expansion of the Weibull PDF and CDF (the CDF comes into play when considering `C₀(t) = F₀(t) - F₁(t) = 1 - (1 - e^{λ t^c}) = e^{λ t^c}`). The resulting density for the **Weibull counting model** is given by Eq. (11) in [1].

# Notes
The computation of the PDF, CDF, etc. for a Weibull counting model requires approximation of an infinite series, which is highly unstable. We therefore make use of `BigInt` and `BigFloat` for all computations under the hood.
Most of the methods related to `WeibullCountModel` are dependent on the series approximation, and will therefore throw a `WeibullCountModel.DivergedError` upon failing to converge. 

# TODO
- [ ] Implement sampling; maybe do so by sampling from `Weibull` distribution rather than directly from `pdf` of `WeibullCountModel`
- [ ] Certain parts can be made parallel. Whether or not it's sufficiently heavy tasks to make it more efficient, we'll see.
- [ ] Rewrite recursive `weibull_alpha` as while-loop, and allow caching to be optional. With current impl, we run into issues when using auto-grad packages.

# References
[1] Adrian, M., Bradlow, E., Fader, P., & McShane, B., Count models based on weibull interarrival times, CoRR, (),  (2013).
