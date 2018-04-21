# Broyden-Fletcher-Goldfarb-Shanno.
function BFGS!(B::Matrix, y::Vector, s::Vector)
    n, m = size(B)
    Bs = B * s
    B[1:n, 1:m] = B - (Bs * Bs') / dot(s, Bs) + (y * y') / dot(s, y)
end

# Systematic rank one.
function SR1!(B::Matrix, y::Vector, s::Vector)
    n, m = size(B)
    yBs = y - B * s
    B[1:n, 1:m] = B + (yBs * yBs') / (yBs' * s)
end
