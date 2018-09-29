rolt = ((x::UInt64, k::Int64) -> (x << k) | ( x>> (64 - k)))

function next(seed::Array{UInt64, 1})
    result_plus = seed[1] + seed[4]    
    new_seed = copy(seed)
    t = seed[2] << 17
    new_seed[3] = xor(new_seed[3], new_seed[1])
    new_seed[4] = xor(new_seed[4], new_seed[2])
    new_seed[2] = xor(new_seed[2], new_seed[3])
    new_seed[1] = xor(new_seed[1], new_seed[4])
    new_seed[3] = xor(new_seed[3], t)
    new_seed[4] = rolt(new_seed[4], 45)
    return result_plus, new_seed
end

function UNIF(seed::Array{UInt64, 1})
    uint, seed = next!(seed)
    return uint()/0xffffffffffffffff, seed
end

function jump(seed::Array{UInt64, 1}) 
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
    return [s0, s1, s2, s3]
end

using Distributions

function get_normal(seed::Array{UInt64, 1})
    r, new_seed = UNIF(seed)
    return quantile(Normal(), UNIF()), new_seed
end

function get_normal(n::Int64, seed::Array{Float64, 1})
    normals = Array{Float64}(undef, n)
    new_seed = [UInt64(0), UInt64(0), UInt64(0), UInt64(0)]
    for i in 1:n
        nor, new_seed = get_normal(seed)
        normals[i] = nor
    end
    return normals, new_seed
end

function get_normal(frst::Int64, scnd::Int64, seed::Array{UInt64, 1})
    normals = Array{Float64}(undef, frst, scnd)
    new_seed = [UInt64(0), UInt64(0), UInt64(0), UInt64(0)]
    for i in 1:frst
        for j in 1:scnd
            nor, new_seed = get_normal(seed)
            normals[i, j] = nor
        end
    end
    return normals, new_seed
end

function get_normal(tup::Tuple{Int64, Int64}, seed::Array{UInt64, 1})
    frst, scnd = tup
    normals = Array{Float64}(undef, frst, scnd)
    new_seed = [UInt64(0), UInt64(0), UInt64(0), UInt64(0)]
    for i in 1:frst
        for j in 1:scnd
            nor, seed = get_normal(seed)
            normals[i, j] = nor
        end
    end
    return normals, new_seed
end
