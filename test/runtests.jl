using DiscreteChoice
using ForwardDiff
using Base.Test # eventually: using Test

# Rosenbrock.
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;
g(x) = ForwardDiff.gradient(f, x);
H(x) = ForwardDiff.hessian(f, x);

# Gradient descent.
@test optimize(f, g, zeros(2), 100000)[1] ≈ ones(2)

# Newton.
@test optimize(f, g, H, zeros(2))[1] ≈ ones(2)

# Quasi-Newton using inverse BFGS.
#@test optimize(f, g, inv_BFGS, zeros(2), true)[1] ≈ ones(2)

# Trust region using Steihaug-Toint's method and exact Hessian.
@test optimize(f, g, H, tcg, zeros(2))[1] ≈ ones(2)

# Trust region using Steihaug-Toint's method and approximative Hessian using BFGS.
@test optimize(f, g, BFGS!, tcg, zeros(2), true)[1] ≈ ones(2)

# Trust region using Steihaug-Toint's method and approximative Hessian using SR1.
@test optimize(f, g, SR1!, tcg, zeros(2), true)[1] ≈ ones(2)