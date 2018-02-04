using DiscreteChoice
using ForwardDiff
using Base.Test

# write your own tests here
@test 1 == 2

rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(rosenbrock, x);

H(x) = ForwardDiff.hessian(rosenbrock, x);

@test newton(f, g, H, zeros(2)) == ones(2)
