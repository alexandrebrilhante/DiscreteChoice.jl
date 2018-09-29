rolt = ((x::UInt64, k::UInt64) -> (x << k) | (x >> (64 - k)))

function set_seed(s::Array{UInt64})
    global seed[:] = s
end

function next()
    s0 = seed[1]
    s1 = seed[2]
    result = s0 + s1
    s1 = xor(s1, s0)
    seed[1] = xor(rolt(s0, UInt64(24)), s1, (s1 << 16))
    seed[2] = rolt(s1, UInt64(37))
    return result
end

function jump()
    JUMP = [0xdf900294d8f554a5, 0x170865df4b3201fc]
    s0::UInt64 = 0
    s1::UInt64 = 0
    for i in 1:2
        for b in 0:63
            if (JUMP[i] & UInt64(1) << b) != 0
                s0 = xor(s0, seed[1])
                s1 = xor(s1, seed[2])
            end
            next()
        end
    end
    global seed = [s0, s1]
end

function UNIF()
    return next()/0xffffffffffffffff
end

function get_normal(m::Int64 = 1)
    return [quantile(Normal(), UNIF()) for i in 1:m]
end
