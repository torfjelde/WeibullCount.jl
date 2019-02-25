using Test

@testset "Van-Wijngaarden" begin include("test_van_wijngaarden.jl") end
@testset "Weibull Count" begin include("test_weibull_count.jl") end
@testset "Distribution" begin include("test_distribution.jl") end
