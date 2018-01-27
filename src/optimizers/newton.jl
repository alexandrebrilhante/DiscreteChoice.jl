function newton(f::Function, g::Function, h::Function,
                x0::Vector, maxiter::Int64 = 500, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    δ *= δ
    H = eye(n)
    dfx = ones(n)
    g(x, dfx)
    while dot(dfx, dfx) > δ && k < maxiter
        g(x, dfx)
        h(x, H)
        x -= H \ dfx
        k += 1
    end
    return x
end
