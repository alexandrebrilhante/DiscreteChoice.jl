abstract type Individuals end

struct LogitMakerIndividuals <: Individuals
    data::Array{Float64, 2} 
    choice::Int64
    n_sim::Int64
end

struct MixedLogitMakerIndividuals <: Individuals
    data::Array{Float64, 2} 
    choice::Int64
    n_sim::Int64
    seed::Array{UInt64, 1}
end
