# DiscreteChoice.jl

[![Build Status](https://travis-ci.org/brilhana/DiscreteChoice.jl.svg?branch=master)](https://travis-ci.org/brilhana/DiscreteChoice.jl)
[![Coverage Status](https://coveralls.io/repos/brilhana/DiscreteChoice.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/brilhana/DiscreteChoice.jl?branch=master)
[![codecov.io](http://codecov.io/github/brilhana/DiscreteChoice.jl/coverage.svg?branch=master)](http://codecov.io/github/brilhana/DiscreteChoice.jl?branch=master)

Tools for estimating discrete choice models in Julia.

# Solvers

- Gradient descent
- Newton
- Quasi-Newton using BFGS or SR1
- Trust region using exact Hessian or BFGS or SR1

# Installation

```julia
Pkg.clone("https://github.com/brilhana/DiscreteChoice.jl.git")
```

# Usage

```julia
using DiscreteChoice
using ForwardDiff

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(f, x);

H(x) = ForwardDiff.hessian(f, x);

optimize(f, g, H, tcg, zeros(2))
```

# Examples
[Logit](https://github.com/brilhana/DiscreteChoice.jl/blob/master/examples/logit.ipynb)
