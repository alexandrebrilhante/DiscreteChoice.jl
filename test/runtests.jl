# TODO: Tests.

using DiscreteChoice, ForwardDiff
using Base.Test

# Rosenbrock.
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;
rosenbrock_g(x) = ForwardDiff.gradient(rosenbrock, x);
rosenbrock_H(x) = ForwardDiff.hessian(rosenbrock, x);

# Gradient descent.
@test optimize(rosenbrock, rosenbrock_g, zeros(2), 100000)[1] ≈ ones(2)

# Newton.
@test optimize(rosenbrock, rosenbrock_g, rosenbrock_H, zeros(2))[1] ≈ ones(2)

# Trust region with Steihaug-Toint's method and exact Hessian.
@test optimize(rosenbrock, rosenbrock_g, rosenbrock_H, tcg, zeros(2))[2] == 26

# Trust region with Steihaug-Toint's method and approximative Hessian using BFGS.
@test optimize(rosenbrock, rosenbrock_g, BFGS!, tcg, zeros(2), true)[1] ≈ ones(2)

# Trust region with Steihaug-Toint's method and approximative Hessian using SR1.
@test optimize(rosenbrock, rosenbrock_g, SR1!, tcg, zeros(2), true)[1] ≈ ones(2)
