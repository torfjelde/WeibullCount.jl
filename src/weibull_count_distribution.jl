using Distributions
import Distributions: cdf, pdf

@doc raw"""
Weibull count model.

**Weibull count model** is a discrete distribution obtained as the *count model with inter-arrival times assumed to be i.i.d. Weibull distributed* [1, Section 3.1].

In words, a Weibull count model models the *number of events* which occur within some time `t`, where the time *between* each of the events are (i.i.d.) Weibull distributed.

The PDF is given by the following formula
```math
P \big( I(t) = i \big) = \sum_{j = 1}^{\infty} \frac{(-1)^{j + i} (\lambda t^c)^j a_j^i}{\Gamma(cj + 1)}, \quad i = 0, 1, 2, \dots
```
where
```math
\alpha_j^0 = \frac{\Gamma(cj + 1)}{\Gamma(j + 1)}, \quad j = 0, 1, 2, \dots
```
and
```math
\alpha_j^{i + 1} = \sum_{m = i}^{j - 1} \alpha_m^i \frac{\Gamma(cj - cm + 1)}{\Gamma(j - m + 1)}, \quad i = 0, 1, 2, \dots, \quad j = i + 1, i + 2, \dots
```
The derivation can be found in [1, Section 3.1].

## Notes
The computation of the PDF, CDF, etc. for a Weibull counting model requires approximation of an infinite series, which is highly unstable. We therefore make use of `BigInt` and `BigFloat` for all computations under the hood.
Most of the methods related to `WeibullCountModel` are dependent on the series approximation, and will therefore throw a `WeibullCountModel.DivergedError` upon failing to converge. 

## References
[1] Adrian, M., Bradlow, E., Fader, P., & McShane, B., Count models based on weibull interarrival times, CoRR, (),  (2013). 

"""
struct WeibullCountModel <: Distribution{Univariate, Discrete}
    λ::BigFloat
    c::BigFloat
    t::BigFloat
end

WeibullCountModel(λ::Number, c::Number, t::Number) = WeibullCountModel(BigFloat(λ), BigFloat(c), BigFloat(1.0))
WeibullCountModel(λ::Number, c::Number) = WeibullCountModel(λ, c, 1.0)

# PDF
pdf(d::WeibullCountModel, x::Array{BigInt}) = begin
    # extract the actual result of the computation
    weibull_count_pdf.(x, d.λ, d.c, d.t)[1]
end
pdf(d::WeibullCountModel, x::AbstractArray) = pdf.(d, BigInt.(x))

pdf(d::WeibullCountModel, x::BigInt) = weibull_count_pdf(x, d.λ, d.c)

# make explicit to prevent ambiguity from `Distributions.jl`
for T ∈ [Int64, Integer, Real, Number]
    eval(quote pdf(d::WeibullCountModel, x::$T) = pdf(d, BigInt(x)) end)
end

# CDF
cdf(d::WeibullCountModel, x::Array{BigInt}) = begin
    # extract the actual result of the computation
    weibull_count_cdf.(x, d.λ, d.c)[1]
end
cdf(d::WeibullCountModel, x::AbstractArray) = cdf(d, BigInt.(x))

cdf(d::WeibullCountModel, x::BigInt) = weibull_count_cdf(x, d.λ, d.c)
# make explicit to prevent ambiguity from `Distributions.jl`
for T ∈ [Int64, Integer, Real, Number]
    eval(quote cdf(d::WeibullCountModel, x::$T) = cdf(d, BigInt(x)) end)
end
