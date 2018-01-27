using Optim

function gradientdescent(f::Function, g::Function, x0::Vector,
                         h::Float64, maxiter::Int64 = 500, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    while dot(dfx, dfx) > δ && k < maxiter
        fsearch(α) = f(x - α * dfx);
        α = Optim.minimizer(optimize(fsearch, 0, h, GoldenSection()))
        dfx = g(x)
        x -= α*dfx
        k += 1
    end
    return x, k
end
