using WeibullCount

begin
    x = 2
    c = 0.36787944117144233
    λ = 20.08553692
    
    d = WeibullCountProcess(λ, c)

    @test pdf(d, x)
    @test pdf(d, [x])

    @test cdf(d, x)
    @test cdf(d, [x])
end
