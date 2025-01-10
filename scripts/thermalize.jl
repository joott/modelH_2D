cd("/home/jkott/perm/modelH_2D/thermalized")

using JLD2
using CodecZlib

include("../src/modelH.jl")

function main()
    @init_state
    # fft_temp = ArrayType{ComplexType}(undef, (L,L,2)) # required as argument for prethermalize(...)
    arrays = make_temp_arrays(state) # required as argument for thermalize(...)
    mass_id = round(m²,digits=3)
    for i in 1:L
        thermalize(state, arrays, m², L^3)
        @show i
        save_state("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(i)_id_$(seed).jld2", state, m², i)
        if i>1
            rm("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(i-1)_id_$(seed).jld2")
        end
    end
    save_state("thermalized_mass_$(mass_id)_L_$(L)_id_$(seed).jld2", state, m²)
    rm("prelim_thermalized_mass_$(mass_id)_L_$(L)_idx_$(L)_id_$(seed).jld2")
end

main()
