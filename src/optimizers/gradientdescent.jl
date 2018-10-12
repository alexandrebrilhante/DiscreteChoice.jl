# TODO: create interface for optimize.

struct GradientDescent <: Solver end

function armijo(f::Function, dfx::Vector, x::Vector, d::Vector;
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

"""
Gradient descent.
"""
function optimize(f::Function, g::Function, x0::Vector{T},
                  kmax::Int64 = 5000; δ::Float64 = 1e-6) where {T<:Real}
    k::Int64 = 0
    δ *= δ
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    while norm(dfx) > δ && k < kmax
        dfx = g(x)
        α = armijo(f, dfx, x, -dfx)
        x -= α * dfx
        k += 1
    end
    return x, k
end
