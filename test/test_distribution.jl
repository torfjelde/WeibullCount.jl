using WeibullCount

begin
    x = 2
    c = 0.36787944117144233
    λ = 20.08553692
    
    d = WeibullCountModel(λ, c)

    @test pdf(d, x) ≥ 0
    @test all(pdf(d, [x]) .≥ 0)

    @test cdf(d, x) ≥ 0
    @test all(cdf(d, [x]) .≥ 0)
end
