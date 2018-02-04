__precompile__(true)

module DiscreteChoice
    export gradientdescent,
           newton,
           trustregion

include("optimizers/gradientdescent.jl")
include("optimizers/newton.jl")
include("optimizers/trustregion.jl")

end
