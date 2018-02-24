function newton(f::Function, g::Function, H::Function,
                x0::Vector, kmax::Int64 = 1000, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    d2fx = eye(n, n)
    while norm(dfx) > δ && k < kmax
        g!(x, dfx)
        H!(x, d2fx)
        x -= d2fx \ dfx
        k += 1
    end
    return x, k
end
