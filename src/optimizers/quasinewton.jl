# Quasi-Newton.

function optimize(f::Function, g::Function, H::Function, x0::Vector{T},
                  approxh::Bool; kmax::Int64 = 5000, δ::Float64 = 1e-6) where {T<:Real}
    k::Int64 = 0
    δ *= δ
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    d2fx = eye(n, n)
    while norm(dfx) > δ && k < kmax
        prev = x
        dfx = g(x)
        x -= d2fx * dfx
        y = g(x) - dfx
        s = x - prev
        d2fx = H(d2fx, y, s)
        k += 1
    end
    return x, k
end
