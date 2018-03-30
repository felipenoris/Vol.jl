
function discrete_returns(series::Array{Float64, 1}) :: Array{Float64, 1}
	n = length(series)
	r = zeros(n-1)
	for i in 2:n
		r[i-1] = series[i]/series[i-1] - 1
	end
	return r
end

function discrete_returns(series::Array{Float64, 2}) :: Array{Float64, 2}
	(n, num_series) = size(series)
	r = zeros(n-1, num_series)
	for c in 1:num_series
		for i in 2:n
			r[i-1, c] = series[i, c]/series[i-1, c] - 1
		end
	end
	return r
end

function log_returns(series::Array{Float64, 1}) :: Array{Float64, 1}
	n = length(series)
	r = zeros(n-1)
	for i in 2:n
		r[i-1] = log(series[i]/series[i-1])
	end
	return r
end

function log_returns(series::Array{Float64, 2}) :: Array{Float64, 2}
	(n, num_series) = size(series)
	r = zeros(n-1, num_series)
	for c in 1:num_series
		for i in 2:n
			r[i-1, c] = log(series[i, c]/series[i-1, c])
		end
	end
	return r
end
