__precompile__(true)

module DiscreteChoice
    export gradientdescent,
           newton,
           trustregion,
           hessian
    end

include("optimizers/gradientdescent.jl")
include("optimizers/newton.jl")
include("optimizers/trustregion.jl")
include("utils/hessian.jl")

end
