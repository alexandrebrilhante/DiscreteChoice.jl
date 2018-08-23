module DiscreteChoice

using LinearAlgebra

import ForwardDiff

export optimize,
       armijo,
       tcg,
       BFGS!,
       inv_BFGS,
       SR1!

include("optimizers/gradientdescent.jl")
include("optimizers/newton.jl")
include("optimizers/quasinewton.jl")
include("optimizers/trustregion.jl")
include("utils/hessian.jl")

end # module
