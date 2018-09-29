module DiscreteChoice

using LinearAlgebra

import Distributions, ForwardDiff

export optimize,
       armijo,
       tcg,
       BFGS!,
       inv_BFGS,
       SR1!

include("batch/batch.jl")

include("individuals/individuals.jl")

include("optimizers/ag.jl")
include("optimizers/gradientdescent.jl")
include("optimizers/newton.jl")
include("optimizers/quasinewton.jl")
include("optimizers/rsag.jl")
include("optimizers/trustregion.jl")

include("utils/hessian.jl")
include("utils/rng.jl")
include("utils/xoshiro.jl")
include("utils/xoshiro256plus.jl")

end # module
