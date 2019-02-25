module WeibullCount

export weibull_count_pdf,
    weibull_count_cdf,
    WeibullCountModel,
    pdf, cdf

using SpecialFunctions
using LRUCache

N_TERMS = 300

function loggamma(x)
    log(gamma(x))
end

# By using parametric types, we can use Float32 or BigFloat, and we will return the same type
"""
    `van_wijngaarden_transformation(a::AbstractArray{T})::T where T <: AbstractFloat`

Computes the alternating sum of the terms `a` (assuming first term to be have positive coefficient).

This method is an example ofa "series acceleration", and will in most cases increase the convergence rate, i.e. requiring fewer terms
to obtain converge for a given tolerance, as compared to just computing the partial sum.
Van-Wijngaarden transformation is a particular implementation of the Euler transformation of the alternating series given by terms `a`.
"""
function van_wijngaarden_transformation(a::AbstractArray{T})::T where T <: AbstractFloat
    n = length(a)
    
    first_row::Array{T} = zeros(n)
    first_row[1] = a[1]
    
    for i = 2:n
        first_row[i] = first_row[i - 1] + (-1)^(i - 1) * a[i]
    end
    
    row = first_row
    
    for i = 1:Int(floor(2 * n / 3))
        row = 0.5 * (first_row[1:end - 1] + first_row[2:end])
    end
    
    return row[end]
end


function _weibull_alpha(j::Integer, i::Integer, c::Number)::Number
    if i == 0
        return exp(loggamma(c * j + 1) - loggamma(j + 1))
    end
    
    return sum(exp(log(weibull_alpha(m, i - 1, c)) + loggamma(c * j - c * m + 1) - loggamma(j - m + 1))
    for m = i - 1:j - 1)
end


const WEIBULL_ALPHA_CACHE_SIZE = 10000
const _weibull_alpha_cache = LRU{Tuple{Integer, Integer, Number}, Number}(WEIBULL_ALPHA_CACHE_SIZE)

cache_clear!() = empty!(_weibull_alpha_cache)
cache_size!(n) = resize!(_weibull_alpha_cache, n)
cache_info() = length(_weibull_alpha_cache)

"""
    `weibull_alpha(j::Integer, i::Integer, c::Number)::Number`

This method used repeadly 
"""
function weibull_alpha(j::Integer, i::Integer, c::Number)::Number
    @get! _weibull_alpha_cache (j, i, c) _weibull_alpha(j, i, c)
end

# HACK: only have cache for the highest
# const _weibull_alpha_cache = LRU{Tuple{BigInt, BigInt, BigFloat}, BigFloat}(32000)
# function weibull_alpha(j::BigInt, i::BigInt, c::BigFloat)::BigFloat
#     @get! _weibull_alpha_cache (j, i, c) _weibull_alpha(j, i, c)
# end

"DEPRECATED."
function weibull_count_pdf_old(x::T1, λ::T2, c::T2, t::T2; k::T1=convert(T1, N_TERMS), transform::Bool=true)::T2 where {T1 <: Integer, T2 <: AbstractFloat}
    js = collect(x:x + k)
    
    terms::Array{T2} = exp.(js .* (log(λ) .+ c * log(t)) .+ log.(weibull_alpha.(js, x, c)) .- loggamma.(c .* js .+ 1))
    
    if transform
        s = van_wijngaarden_transformation(terms)
        
        # TODO: start with fewer terms, and then keep adding until we have "convergence"
    else
        s = sum((-1).^(x .+ js) .* terms)
    end
    
    return s
end

"""
    `weibull_count_pdf_approx(x::T1, λ::T2, c::T2, t::T2; k::T1 = convert(T1, $N_TERMS), transform::Bool = true, increment::T1 = convert(T1, 10), tol::T2 = convert(T2, 0.001))::Tuple{T2, Bool} where {T1 <: Integer, T2 <: AbstractFloat}`

Computes approximation of Eq. 11 in [1] using partial sum and, if `transform = true`, Van-Wijngaarden transformation.

## References
[1] Adrian, M., Bradlow, E., Fader, P., & McShane, B., Count models based on weibull interarrival times, CoRR, (),  (2013). 

"""
function weibull_count_pdf_approx(x::T1, λ::T2, c::T2, t::T2; k::T1 = convert(T1, N_TERMS), transform::Bool = true, increment::T1 = convert(T1, 10), tol::T2 = convert(T2, 0.001))::Tuple{T2, Bool} where {T1 <: Integer, T2 <: AbstractFloat}
    # endpoint
    maximum = x + k
    
    js = collect(x:x + increment)
    
    terms::Array{T2} = exp.(js .* (log(λ) .+ c * log(t)) .+ log.(weibull_alpha.(js, x, c)) .- loggamma.(c .* js .+ 1))

    if transform
        s = van_wijngaarden_transformation(terms)
    else
        s = sum((-1).^(x .+ js) .* terms)
    end
        
    converged = false
    valid = false

    # incremently add more terms
    for i = increment:increment:maximum
        js = collect(x + i + 1: min(maximum, x + i + increment))
        append!(terms, exp.(js .* (log(λ) .+ c * log(t)) .+ log.(weibull_alpha.(js, x, c)) .- loggamma.(c .* js .+ 1)))

        # print(i)
        # print(js)

        # compute difference
        if transform
            s_new = van_wijngaarden_transformation(terms)
        else
            s_new = sum((-1).^(x .+ js) .* terms)
        end
        
        Δ = abs(s_new - s)
        
        # println("Delta: $Δ")
        # println("s: $s")

        converged = Δ ≤ tol
        valid = (s ≥ T2(0.0)) & (s ≤ T2(1.0))

        if converged & valid
            break
        end

        # update
        s = s_new
    end

    # if !converged
    #     println("failed to converge")
    # end
    # if !valid
    #     println("invalid probability")
    # end

    return s, (valid & converged)
end

struct DivergedError{T1, T2} <: Exception where {T1 <: Integer, T2 <: AbstractFloat}
    x::T1
    λ::T2
    c::T2
    t::T2
    k::T1
    transform::Bool
    increment::T1
    tol::T2
    result::T2
end

Base.showerror(io::IO, e::DivergedError) = print(io, "series approximation failed to converge; try increasing `k`:\nx = $(e.x), λ = $(e.λ), c = $(e.c), t = $(e.t), k = $(e.k), transform = $(e.transform), increment = $(e.increment), tol = $(e.tol), result = $(e.result)")

"""
    `weibull_count_pdf(x::T1, λ::T2, c::T2, t::T2; k::T1 = convert(T1, $N_TERMS), transform::Bool = true, increment::T1 = convert(T1, 10), tol::T2 = convert(T2, 0.001))::T2 where {T1 <: Integer, T2 <: AbstractFloat}`

Computes an approximation to `P(X(t) = x)` for time `t`, for a Weibull counting model with parameters `λ` and `c`.

The main computation is performed by `WeibullCount.weibull_count_pdf_approx`

## See also
`WeibullCount.weibull_count_pdf_approx`

"""
function weibull_count_pdf(x::T1, λ::T2, c::T2, t::T2; k::T1 = convert(T1, N_TERMS), transform::Bool = true, increment::T1 = convert(T1, 10), tol::T2 = convert(T2, 0.001))::T2 where {T1 <: Integer, T2 <: AbstractFloat}
    p, converged = weibull_count_pdf_approx(x, λ, c, t; k = k, transform = transform, increment = increment, tol = tol)

    if !converged
        throw(DivergedError(x, λ, c, t, k, transform, increment, tol, p))
    end

    return p
end

function weibull_count_pdf(x::T1, λ::T2, c::T2; k::T1 = convert(T1, N_TERMS), transform::Bool = true, increment::T1 = convert(T1, 10), tol::T2 = convert(T2, 0.001))::T2 where {T1 <: Integer, T2 <: AbstractFloat}
    weibull_count_pdf(x, λ, c, T2(1.0); k = k, transform = transform, increment = increment, tol = tol)
end

function weibull_count_cdf(x::T1, λ::T2, c::T2; k::T1=convert(T1, N_TERMS))::T2 where {T1 <: Integer, T2 <: AbstractFloat}
    if x < 0
        return 0
    elseif x == 0
        return sum(weibull_count_pdf(0, λ, c; k=k))
    end
    
    return sum(weibull_count_pdf(y, λ, c; k=k) for y = 0:x)
end


# includes
include("weibull_count_distribution.jl")

end # module
