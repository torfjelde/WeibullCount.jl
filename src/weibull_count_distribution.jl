using Distributions
import Distributions: cdf, pdf

struct WeibullCountProcess <: Distribution{Univariate, Discrete}
    λ::BigFloat
    c::BigFloat
    t::BigFloat
end

WeibullCountProcess(λ::Number, c::Number, t::Number) = WeibullCountProcess(BigFloat(λ), BigFloat(c), BigFloat(1.0))
WeibullCountProcess(λ::Number, c::Number) = WeibullCountProcess(λ, c, 1.0)

# PDF
pdf(d::WeibullCountProcess, x::Array{BigInt}) = begin
    # extract the actual result of the computation
    weibull_count_pdf.(x, d.λ, d.c, d.t)[1]
end
pdf(d::WeibullCountProcess, x::AbstractArray) = pdf.(d, BigInt.(x))

pdf(d::WeibullCountProcess, x::BigInt) = weibull_count_pdf(x, d.λ, d.c)

# make explicit to prevent ambiguity from `Distributions.jl`
for T ∈ [Int64, Integer, Real, Number]
    eval(quote pdf(d::WeibullCountProcess, x::$T) = pdf(d, BigInt(x)) end)
end

# CDF
cdf(d::WeibullCountProcess, x::Array{BigInt}) = begin
    # extract the actual result of the computation
    weibull_count_cdf.(x, d.λ, d.c)[1]
end
cdf(d::WeibullCountProcess, x::AbstractArray) = cdf(d, BigInt.(x))

cdf(d::WeibullCountProcess, x::BigInt) = weibull_count_cdf(x, d.λ, d.c)
# make explicit to prevent ambiguity from `Distributions.jl`
for T ∈ [Int64, Integer, Real, Number]
    eval(quote cdf(d::WeibullCountProcess, x::$T) = cdf(d, BigInt(x)) end)
end
