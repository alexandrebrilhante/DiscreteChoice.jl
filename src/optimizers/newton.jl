# Newton.

function optimize(f::Function, g::Function, H::Function,
                  x0::Vector{T}; kmax::Int64 = 5000, δ::Float64 = 1e-6) where {T<:Real}
    k::Int64 = 0
    δ *= δ
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    d2fx = eye(n, n)
    while norm(dfx) > δ && k < kmax
        dfx = g(x)
        d2fx = H(x)
        x -= d2fx \ dfx
        k += 1
    end
    return x, k
end