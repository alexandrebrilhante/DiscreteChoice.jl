using Optim

function gradientdescent(f::Function, g::Function, x0::Vector{Float64},
                         h::Float64, maxiter::Int64 = 1000, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    fsearch(α::Float64) = f(x - α * dfx);
    while norm(dfx) > δ && k < maxiter
        α = Optim.minimizer(optimize(fsearch, 0.0, h, GoldenSection()))
        dfx = g(x)
        x -= α * dfx
        k += 1
    end
    return x, k
end
