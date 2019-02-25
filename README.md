# Weibull count model.

**Weibull count model** is a discrete distribution obtained as the *count model with inter-arrival times assumed to be i.i.d. Weibull distributed* [1, Section 3.1].

That is, a Weibull count model models the *number of events* which occur within some time `t`, where the time *between* each of the events is (i.i.d.) Weibull distributed.

## 

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

    F₁(t) = ∫ f(s) d(s),    s ∈ [0, t]
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

# References
[1] Adrian, M., Bradlow, E., Fader, P., & McShane, B., Count models based on weibull interarrival times, CoRR, (),  (2013).
