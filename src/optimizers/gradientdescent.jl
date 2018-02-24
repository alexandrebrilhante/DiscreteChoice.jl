function armijo(f::Function, dfx::Vector, x::Vector, d::Vector,
                α::Float64 = 1.0, β1::Float64 = 1e-4, β2::Float64 = 0.9)
    s = β1 * dot(dfx, d)
    fx = f(x)
    fxcand = f(x + α * d)
    while fxcand > fx + α * s
        α *= β2
        fxcand = f(x + α * d)
    end
    return α
end

function gradientdescent(f::Function, g::Function, x0::Vector,
                         h::Float64, kmax::Int64 = 1000, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    while norm(dfx) > δ && k < kmax
        g!(x, dfx)
        α = armijo(f, dfx, x, -dfx)
        x -= α * dfx
        k += 1
    end
    return x, k
end
