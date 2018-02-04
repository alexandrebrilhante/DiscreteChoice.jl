function newton(f::Function, g::Function, H::Function,
                x0::Vector{Float64}, maxiter::Int64 = 1000, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    hessian = eye(n)
    while norm(dfx) > δ && k < maxiter
        dfx = gradient = g(x)
        hessian = H(x)
        x -= hessian \ gradient
        k += 1
    end
    return x
end
