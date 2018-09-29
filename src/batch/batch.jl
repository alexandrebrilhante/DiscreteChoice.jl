abstract type Batch end

function iterate(b::Batch, state::Tuple{Int64, Int64} = (1, 1))
    if b.done > b.n
        b.done = 1
        return nothing
    else
        b.done += 1
        try
            x = state[1]
            y = state[2]
            n_alt = b.n_alt(x)
            n_sim = b.n_sim(x)
            data_i = b.data[y:y + n_alt - 1, :]
            new_state = (x + 1, y + n_alt)
            return LogitMakerIndividuals(data_i, 1, n_sim), new_state
        catch
            x = 1
            y = 1
            n_alt = b.n_alt(x)
            n_sim = b.n_sim(x)
            data_i = b.data[y:y + n_alt - 1, :]
            new_state = (x + 1, y + n_alt)
            return LogitMakerIndividuals(data_i, 1, n_sim), new_state
        end
    end
end
