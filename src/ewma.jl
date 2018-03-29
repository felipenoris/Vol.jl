
"""
    ewma(series, λ) -> Float64

Canonical definition of Exponentially Weighted Moving Average.
"""
function ewma_canonical(series::Vector{Float64}, λ::Float64) :: Float64
    n = length(series)

    local sumprod::Float64 = 0.0
    local sumweights::Float64 = 0.0

    for i in 1:n
        weight = λ^(i-1)
        sumprod += series[i] * weight
        sumweights += weight
    end
    return sumprod / sumweights
end

"""
Same as ewma_canonical, but squares returns and returns sqrt of the result.
"""
function ewma(returns::Vector{Float64}, λ::Float64) :: Float64
    n = length(returns)

    local sumprod::Float64 = 0.0
    local sumweights::Float64 = 0.0

    for i in 1:n
        weight = λ^(i-1)
        sumprod += (returns[i]^2) * weight
        sumweights += weight
    end
    return sqrt(sumprod / sumweights)
end
