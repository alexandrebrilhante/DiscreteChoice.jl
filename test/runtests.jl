# TODO: Tests.

using DiscreteChoice, ForwardDiff
using Base.Test

f(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(f, x);

H(x) = ForwardDiff.hessian(f, x);

optimize(f, g, zeros(2))

optimize(f, g, H, zeros(2))

optimize(f, g, H, tcg, zeros(2))

optimize(f, g, BFGS!, tcg, zeros(2), true)

optimize(f, g, SR1!, tcg, zeros(2), true)
