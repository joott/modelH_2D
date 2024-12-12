cd("/home/jkott/perm/modelH_2D/thermalized")

using JLD2
using CodecZlib

include("src/modelH.jl")

function main()
    @init_state
    arrays = make_temp_arrays(state)
    file = jldopen(init_arg, "r")
    idx0 = file["i"]

    mass_id = round(file["m²"],digits=3)
    for i in idx0+1:L
        thermalize(state, arrays, mass_id, L^3)
        @show i
        save_state("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(i)_id_$(seed).jld2", state, mass_id, i)
        if i>idx0+1
            rm("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(i-1)_id_$(seed).jld2")
        else
            rm(init_arg)
        end
    end
    save_state("thermalized_mass_$(mass_id)_L_$(L)_id_$(seed).jld2", state, mass_id)
    rm("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(L)_id_$(seed).jld2")
end

main()
