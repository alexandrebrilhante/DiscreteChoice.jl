# DiscreteChoice.jl

[![Build Status](https://travis-ci.org/brilhana/DiscreteChoice.jl.svg?branch=master)](https://travis-ci.org/brilhana/DiscreteChoice.jl)
[![Coverage Status](https://coveralls.io/repos/brilhana/DiscreteChoice.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/brilhana/DiscreteChoice.jl?branch=master)
[![codecov.io](http://codecov.io/github/brilhana/DiscreteChoice.jl/coverage.svg?branch=master)](http://codecov.io/github/brilhana/DiscreteChoice.jl?branch=master)

Discrete choice modeling and estimation using stochastic approximation methods in Julia.

## Models

* `Logit`
* `MixedLogit`

## Solvers

* `AG` (Accelerated Gradient)
* `GradientDescent`
* `Newton`
* `QuasiNewton`
* `RSAG` (Randomized Stochastic Accelerated Gradient)
* `SteihaugToint`

## Installation

```julia
(v1.0) pkg> add https://github.com/brilhana/DiscreteChoice.jl
```

## Examples
* [Logit](examples/logit.ipynb)
