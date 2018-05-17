__precompile__(true)

module DiscreteChoice

    using ForwardDiff

    export optimize,
           armijo,
           tcg,
           BFGS!,
           SR1!

    include("optimizers/gradientdescent.jl")
    include("optimizers/newton.jl")
    include("optimizers/quasinewton.jl")
    include("optimizers/trustregion.jl")
    include("utils/hessian.jl")

end
