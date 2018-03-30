
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
