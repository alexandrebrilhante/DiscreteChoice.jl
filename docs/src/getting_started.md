# Getting Started

`DiscreteChoice` is a unregistered package for the moment. To add it to your Julia packages, simply do the following in REPL:

```julia
(v1.0) pkg> add https://github.com/brilhana/DiscreteChoice.jl
```

The function, and its gradient and Hessian (depending on the method) must be supplied to the solver by the user. An easy way to do it is to use `ForwardDiff`.

```julia
using DiscreteChoice
using ForwardDiff

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(f, x);

H(x) = ForwardDiff.hessian(f, x);

# Trust region using Steihaug-Toint.
optimize(f, g, H, tcg, zeros(2))
```
