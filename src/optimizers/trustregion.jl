# Trust region method using the truncated conjugate gradient.

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
    k::Int64
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
    end
    return false
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


# Steihaug-Toint's method.
function tcg(H::Matrix, g::Vector, Δ::Float64)
    n::Int64 = length(g)
    s = zeros(n)
    normg0 = norm(g)
    v = g
    d = -v
    norm2d = gv = dot(g, v)
    norm2s::Float64 = 0.0
    sMd::Float64 = 0.0
    k::Int64 = 0
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
    end
    return false
end

function optimize(f::Function, g::Function, H::Function, step::Function,
                  x0::Vector{T}, approxh::Bool = false, kmax::Int64 = 5000) where {T<:Real}
    tr = trdefaults()
    state::TrustRegionState = TrustRegionState()
    state.k = 0
    state.x = x0
    state.g = g(state.x)
    state.Δ = 1.0
    state.δ *= state.δ
    n = length(state.x)
    d2fx = eye(n, n)
    fx = f(state.x)
    if approxh
        y = zeros(n)
        gcand = zeros(n)
    else
        d2fx = H(state.x)
    end

    function model(s::Vector, g::Vector, H::Matrix)
        return dot(s, g) + 0.5 * dot(s, H * s)
    end

    while norm(state.g) > state.δ && state.k < kmax
        state.step = step(d2fx, state.g, state.Δ)
        state.xcand = state.x + state.step
        fcand = f(state.xcand)
        state.ρ = (fcand - fx) / model(state.step, state.g, d2fx)
        if approxh
            gcand = g(state.xcand)
            y = gcand - state.g
            d2fx = H(d2fx, y, state.step)
        end
        if acceptcandidate!(state, tr)
            state.x = copy(state.xcand)
            if !approxh
                state.g = g(state.x)
                d2fx = H(state.x)
            else
                state.g = copy(gcand)
            end
            fx = fcand
        end
        updateradius!(state, tr)
        state.k += 1
    end
    return state.x, state.k
end
