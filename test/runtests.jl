
using Base.Test
import Vol

return_series = 0.01 * collect(1:50)
lambda = 0.94

@test isapprox(Vol.ewma(return_series, lambda), sqrt(Vol.ewma_canonical(return_series.^2, lambda)))
@test isapprox(Vol.ewma(0.01 * ones(100), lambda), 0.01)

r = [ 1. 2. 3.;
      3. 4. 5.;
      5. 6. 7.;
      10. 20. 30.]

@test issymmetric(Vol.ewma(r, lambda))

single_serie = [1.0, 3.0, 10.0]
single_discrete_returns = zeros(2)
single_discrete_returns[1] = (3.0-1.0)/1.0
single_discrete_returns[2] = (10.0-3.0)/3.0
@test isapprox(Vol.discrete_returns(single_serie), single_discrete_returns)

single_log_returns = zeros(2)
single_log_returns[1] = log(3.0) - log(1.0)
single_log_returns[2] = log(10.0) - log(3.0)
@test isapprox(Vol.log_returns(single_serie), single_log_returns)

two_series = [ 1.0 10.0;
               3.0 5.0;
               2.0 15.0]

two_discrete_returns=zeros(2,2)
two_discrete_returns[1,1] = (3.0-1.0)/1.0
two_discrete_returns[2,1] = (2.0-3.0)/3.0
two_discrete_returns[1,2] = (5.0-10.0)/10.0
two_discrete_returns[2,2] = (15.0-5.0)/5.0
@test isapprox(Vol.discrete_returns(two_series), two_discrete_returns)

two_log_returns=zeros(2,2)
two_log_returns[1,1] = log(3.0) - log(1.0)
two_log_returns[2,1] = log(2.0) - log(3.0)
two_log_returns[1,2] = log(5.0) - log(10.0)
two_log_returns[2,2] = log(15.0) - log(5.0)
@test isapprox(Vol.log_returns(two_series), two_log_returns)

@test isapprox(Vol.discrete_returns(r)[1:end, 1], Vol.discrete_returns(r[1:end, 1]))
@test isapprox(Vol.log_returns(r)[1:end, 1], Vol.log_returns(r[1:end, 1]))

single_return = [0.05, 0.01, 0.10, -0.02, -0.10, 0.0]
lambda = 0.94
ewma_canonical_single_return = (0.0 * 1 + (-0.10)*lambda + (-0.02)*(lambda^2) + 0.10*(lambda^3) + 0.01*(lambda^4) + 0.05*(lambda^5)) / (1 + lambda + lambda^2 + lambda^3 + lambda^4 + lambda^5)
@test isapprox(Vol.ewma_canonical(single_return, lambda), ewma_canonical_single_return)

ewma_single_return = sqrt( ((0.0^2) * 1 + ((-0.10)^2)*lambda + ((-0.02)^2)*(lambda^2) + (0.10^2)*(lambda^3) + (0.01^2)*(lambda^4) + (0.05^2)*(lambda^5)) / (1 + lambda + lambda^2 + lambda^3 + lambda^4 + lambda^5))
@test isapprox(Vol.ewma(single_return, lambda), ewma_single_return)
