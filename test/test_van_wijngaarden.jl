using WeibullCount

begin
    """
    Tests the Van Wijngaarden transformation on the Leibniz formula for π:

        1 - (1 / 3) + (1 / 5) - (1 / 7) + ... = π / 4

    """
    
    a = 1 ./ (1:2:100)
    coeffs = ones(length(a))
    coeffs[2:2:length(a)] .= -1  # alternating minus-sign
    
    true_sum = π / 4
    partial = sum(a .* coeffs)
    trans = WeibullCount.van_wijngaarden_transformation(a)
    
    @test abs(true_sum - partial) .> abs(true_sum - trans)
end
