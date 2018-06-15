# Broyden-Fletcher-Goldfarb-Shanno.
function BFGS!(B::Matrix, y::Vector, s::Vector)
    n, m = size(B)
    Bs = B * s
    B[1:n, 1:m] = B - (Bs * Bs') / dot(s, Bs) + (y * y') / dot(s, y)
end

# Inverse Broyden-Fletcher-Goldfarb-Shanno.
function inv_BFGS(B::Matrix, y::Vector, s::Vector)
    ys = dot(y, s)
    By = B * y
    return B + 1.0 / ys * ((ys + dot(y, By)) / ys * (s * s') - By * s' - s * By')
end

# Systematic rank one.
function SR1!(B::Matrix, y::Vector, s::Vector)
    n, m = size(B)
    yBs = y - B * s
    B[1:n, 1:m] = B + (yBs * yBs') / (yBs' * s)
end
