struct AcceleratedGradient <: Solver end

"""
Accelerated gradient.
"""
function optimize(f::Function, g::Function, x0::Vector{T}, m::AcceleratedGradient, kmax::Int64 = 5000;
                  δ::Float64 = 1e-6) where {T<:Real}
    k::Int64 = 0
    δ *= δ
    x = x0
    y = x0
    n::Int64 = length(x)
    dfx = ones(n)
    while norm(dfx) > δ && k < kmax
        x_prev = x
        dfx = g(x)
        y = x - dfx
        x = y + (k / (k + 3)) * (y - y_prev)
        k += 1
    end
    return x, k
end
