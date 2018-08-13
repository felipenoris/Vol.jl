
"""
    ewma_canonical(series, λ) -> Float64

Canonical definition of Exponentially Weighted Moving Average.
"""
function ewma_canonical(series::Vector{Float64}, λ::Float64) :: Float64
    n = length(series)

    local sumprod::Float64 = 0.0
    local sumweights::Float64 = 0.0

    for (i, serie_element) in enumerate(series)
        weight = λ^(n-i)
        sumprod += serie_element * weight
        sumweights += weight
    end

    return sumprod / sumweights
end

"""
    ewma(returns, λ) -> Float64

Same as `ewma_canonical`, but squares returns and applies sqrt on the result.
Use this function to estimate the volatility of a return series with the EWMA model.
"""
function ewma(returns::Vector{Float64}, λ::Float64) :: Float64
    n = length(returns)

    local sumprod::Float64 = 0.0
    local sumweights::Float64 = 0.0

    for (i, return_element) in enumerate(returns)
        weight = λ^(n-i)
        sumprod += (return_element^2) * weight
        sumweights += weight
    end

    return sqrt(sumprod / sumweights)
end

"""
    ewma_cov(returns::Array{Float64, 2}, λ::Float64) -> Array{Float64, 2}

In the case of a matrix of returns, considers that each column is a time series of returns.
The result is the covariance matrix.
"""
function ewma_cov(returns::Array{Float64, 2}, λ::Float64) :: Array{Float64, 2}
    (rows, cols) = size(returns)
    covmatrix = Array{Float64, 2}(undef, cols, cols)

    @inbounds for j in 1:cols, i in 1:j
        local sumprod::Float64 = 0.0
        local sumweights::Float64 = 0.0

        for r in 1:rows
            weight = λ^(rows-r)
            sumprod += weight * returns[r, i] * returns[r, j]
            sumweights += weight
        end

        covmatrix[i,j] = sumprod / sumweights

        if i != j
            covmatrix[j,i] = covmatrix[i,j]
        end
    end

    return covmatrix
end

"""
    ewma(returns::Array{Float64, 2}, λ::Float64) -> Array{Float64, 2}

In the case of a matrix of returns, considers that each column is a time series of returns.
The result is the correlation matrix.
"""
function ewma(returns::Array{Float64, 2}, λ::Float64) :: Array{Float64, 2}
    (rows, cols) = size(returns)
    covmatrix = ewma_cov(returns, λ)
    @assert size(covmatrix) == (cols, cols)

    # vols is the Volatility vector
    # We can get them from the sqrt of the diagonal of the cov matrix
    vols = sqrt.(LinearAlgebra.diag(covmatrix))

    # turns covmatrix into a correlation matrix
    @inbounds for j in 1:cols, i in 1:j
        covmatrix[i,j] = covmatrix[i,j] / (vols[i] * vols[j])

        if i != j
            covmatrix[j,i] = covmatrix[i,j]
        end
    end

    return covmatrix
end
