function newton(f::Function, g::Function, h::Function,
                x0::Vector, maxiter::Int64, tol::Float64 = 1e-6)
    k = 0
    x = x0
    n = length(x)
    tol *= tol
    H = eye(n)
    dfx = ones(n)
    g(x, dfx)
    while dot(dfx, dfx) > tol && k < maxiter
        k += 1
        g(x, dfx)
        h(x, H)
        x -= H \ dfx
    end
    return x
end
