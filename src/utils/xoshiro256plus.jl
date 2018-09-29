rolt = ((x::UInt64, k::Int64) -> (x << k) | (x >> (64 - k)))

function set_seed(s::Array{UInt64})
    global seed[:] = s
end

function next()
    result_plus = seed[1] + seed[4]    
    t = seed[2] << 17
    seed[3] = xor(seed[3], seed[1])
    seed[4] = xor(seed[4], seed[2])
    seed[2] = xor(seed[2], seed[3])
    seed[1] = xor(seed[1], seed[4])
    seed[3] = xor(seed[3], t)
    seed[4] = rolt(seed[4], 45)
    return result_plus
end

function jump()
    JUMP = [0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c]
    s0::UInt64 = 0
    s1::UInt64 = 0
    s2::UInt64 = 0
    s3::UInt64 = 0
    for i in 1:4
        for b in 0:63
            if (JUMP[i] & (UInt64(1) << b)) != 0
                s0 = xor(s0, seed[1])
                s1 = xor(s1, seed[2])
                s3 = xor(s2, seed[3])
                s4 = xor(s3, seed[4])
            end
            next()
        end
    end
    set_seed([s0, s1, s2, s3]) 
end

function long_jump()
    LONG_JUMP = [0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635]
    s0 = UInt64(0)
    s1 = UInt64(0)
    s2 = UInt64(0)
    s3 = UInt64(0)
    for i in 1:4
        for b in 0:63
            if (LONG_JUMP[i] & (UInt64(1) << b)) != 0
                s0 = xor(s0, seed[1])
                s1 = xor(s1, seed[2])
                s3 = xor(s2, seed[3])
                s4 = xor(s3, seed[4])
            end
            next()
        end
    end
    set_seed([s0, s1, s2, s3]) 
end

function UNIF()
    return next()/0xffffffffffffffff
end

function get_normal()
    quantile(Normal(), UNIF())
end

function get_normal(n::Int64)
    return [quantile(Normal(), UNIF()) for i in 1:n]
end

function get_normal(frst::Int64, scnd::Int64)
    return [quantile(Normal(), UNIF()) for i in 1:frst, j in 1:scnd]
end

function get_normal(tup::Tuple{Int64, Int64})
    return [quantile(Normal(), UNIF()) for i in 1:tup[1], j in 1:tup[2]]
end
