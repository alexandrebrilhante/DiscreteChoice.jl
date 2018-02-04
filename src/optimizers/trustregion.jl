struct TrustRegion{T<:Real}
    η1::T
    η2::T
    γ1::T
    γ2::T
end

function trdefaults()
    return TrustRegion(0.01, 0.9, 0.5, 0.5)
end

mutable struct TrustRegionState
    iter::Int64
    x::Vector{Float64}
    xcand::Vector{Float64}
    g::Vector{Float64}
    step::Vector{Float64}
    Δ::Float64
    ρ::Float64
    δ::Float64

    function TrustRegionState()
        state = new()
        state.δ = 1e-6
        return state
    end
end

function acceptcandidate!(state::TrustRegionState, tr::TrustRegion)
    if state.ρ >= tr.η1
        return true
    else
        return false
    end
end

function updateradius!(state::TrustRegionState, tr::TrustRegion)
    if state.ρ >= tr.η2
        stepnorm = norm(state.step)
        state.Δ = min(1e20, max(4 * stepnorm, state.Δ))
    elseif state.ρ >= tr.η1
        state.Δ *= tr.γ2
    else
        state.Δ *= tr.γ1
    end
end

function cg(A::Matrix{Float64}, b::Vector{Float64},
            x0::Vector{Float64}, δ::Float64 = 1e-6)
    n = length(x0)
    x = x0
    g = b + A * x
    d = -g
    k::Int64 = 0
    δ *= δ
    while dot(g, g) > δ
        Ad = A * d
        normd = dot(d, Ad)
        α = -dot(d, g) / normd
        x += α * d
        g = b + A *x
        γ = dot(g, Ad) / normd
        d = -g + γ * d
        k += 1
    end
    normd = dot(d, A * d)
    α = -dot(d, g) / normd
    x += α * d
    return x
end

function truncatedcg(g::Vector{Float64}, H::Matrix{Float64}, Δ::Float64)
    n::Int64 = length(g)
    s = zeros(n)
    normg0 = norm(g)
    v = g
    d = -v
    gv = dot(g, v)
    norm2d = gv
    norm2s::Float64 = 0.0
    sMd::Float64 = 0.0
    k::Int64 = 0
    Δ *= Δ
    while !stopcg(norm(g), normg0, k, n)
        Hd = H * d
        κ = dot(d, Hd)
        if κ <= 0
            σ = (-sMd + sqrt(sMd * sMd + norm2d * (Δ - dot(s, s)))) / norm2d
            s += σ * d
            break
        end
        α = gv / κ
        norm2s += α * (2 * sMd + α * norm2d)
        if norm2s >= Δ
            σ = (-sMd + sqrt(sMd * sMd + norm2d * (Δ - dot(s, s)))) / norm2d
            s += σ * d
            break
        end
        s += α * d
        g += α * Hd
        v = g
        newgv = dot(g, v)
        β = newgv / gv
        gv = newgv
        d = -v + β * d
        sMd = β * (sMd + α * norm2d)
        norm2d = gv + β * β * norm2d
        k += 1
    end
    return s
end

function stopcg(normg::Float64, normg0::Float64, k::Int64, kmax::Int64)
    χ::Float64 = 0.1
    θ::Float64 = 0.5
    if k == kmax || normg <= normg0 * min(χ, normg0^θ)
        return true
    else
        return false
    end
end

function tr(f::Function, g::Function, H::Function,
            step::Function, x0::Vector{Float64}, maxiter::Int64 = 1000,
            approxH::Bool = false)
    state::TrustRegionState = TrustRegionState()
    tr = trdefaults()
    state.iter = 0
    state.x = x0
    n::Int64 = length(x0)
    δ = state.δ * state.δ
    state.g = zeros(n)
    H = eye(n, n)
    fx = f(x0)
    g(x0, state.g)
    state.Δ = 0.1 * norm(state.g)
    if approxH
        y = zeros(n)
        gcand = zeros(n)
    else
        H(x0, H)
    end

    function model(s::Vector{Float64}, g::Vector{Float64}, H::Matrix{Float64})
        return dot(s, g) + 0.5 * dot(s, H * s)
    end

    while dot(state.g, state.g) > δ && state.iter < maxiter
        if step == truncatedcg
            state.step = step(state.g, H, state.Δ)
        elseif step == cg
            state.step = step(H, state.g, x0)
        end
        state.xcand = state.x + state.step
        fcand = f(state.xcand)
        state.ρ = (fcand - fx) / model(state.step, state.g, H)
        if approxH
            g(state.xcand, gcand)
            y = gcand - state.g
            H = H(H, y, state.step)
        end
        if acceptcandidate!(state, tr)
            state.x = copy(state.xcand)
            if !approxH
                g(state.x, state.g)
                H(state.x, H)
            else
                state.g = copy(gcand)
            end
            fx = fcand
        end
        updateradius!(state, tr)
        state.iter += 1
    end
    return state.x
end
