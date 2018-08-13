
function discrete_returns(series::Vector{Float64}) :: Vector{Float64}
    n = length(series)
    r = Vector{Float64}(undef, n-1)

    @inbounds for i in 2:n
        r[i-1] = series[i]/series[i-1] - 1
    end

    return r
end

function discrete_returns(series::Array{Float64, 2}) :: Array{Float64, 2}
    (n, num_series) = size(series)
    r = Array{Float64, 2}(undef, n-1, num_series)

    @inbounds for c in 1:num_series, i in 2:n
        r[i-1, c] = series[i, c]/series[i-1, c] - 1
    end

    return r
end

function log_returns(series::Vector{Float64}) :: Vector{Float64}
    n = length(series)
    r = Vector{Float64}(undef, n-1)

    @inbounds for i in 2:n
        r[i-1] = log(series[i]/series[i-1])
    end

    return r
end

function log_returns(series::Array{Float64, 2}) :: Array{Float64, 2}
    (n, num_series) = size(series)
    r = Array{Float64, 2}(undef, n-1, num_series)

    @inbounds for c in 1:num_series, i in 2:n
        r[i-1, c] = log(series[i, c]/series[i-1, c])
    end

    return r
end
