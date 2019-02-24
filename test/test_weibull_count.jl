using WeibullCount


function test_weibull_count_pdf(x, λ, c; k = 200, transformed=true)
    for T1 ∈ [Int64, Int128, BigInt]
        for T2 ∈ [Float64, BigFloat]
            WeibullCount.cache_clear!()
            
            print("$T1 and $T2: ")
            println(@time weibull_count_pdf(T1(x), T2(λ), T2(c), T2(1.0); k = T1(k), transform=transformed))
        end
    end
end

# Ensure `weibull_count_pdf_approx` works with different types.
@testset "Parametric types" begin
    x = 2
    c = 0.36787944117144233
    λ = 20.08553692

    for T1 ∈ [Int64, Int128, BigInt]
        for T2 ∈ [Float64, BigFloat]
            WeibullCount.cache_clear!()
            
            # print("$T1 and $T2: ")
            res, converged = WeibullCount.weibull_count_pdf_approx(T1(x), T2(λ), T2(c), T2(1.0))
            @test isa(res, T2)
        end
    end
end

# Failure of convergence
@testset "DivergedError" begin
    x = 2
    c = 0.36787944117144233
    λ = 20.08553692

    WeibullCount.cache_clear!()
    @test_throws WeibullCount.DivergedError WeibullCount.weibull_count_pdf(x, λ, c)

    WeibullCount.cache_clear!()
    @test WeibullCount.weibull_count_pdf(BigInt(x), BigFloat(λ), BigFloat(c)) ≥ 0
end

# Difficult parameters
@testset "Run-times" begin
    x = 2
    c = 0.36787944117144233
    λ = 20.08553692

    WeibullCount.cache_clear!()
    t1 = @elapsed WeibullCount.weibull_count_pdf_old(BigInt(x), BigFloat(λ), BigFloat(c), BigFloat(1.0); k = BigInt(200))

    WeibullCount.cache_clear!()
    t2 = @elapsed weibull_count_pdf(BigInt(x), BigFloat(λ), BigFloat(c), BigFloat(1.0); k = BigInt(200), increment = BigInt(10), tol = BigFloat(0.001))

    @test t2 < t1
end
