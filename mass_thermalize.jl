cd(@__DIR__)

using JLD2
using CodecZlib

include("src/modelH.jl")

function main()
    @init_state
    fft_temp = ArrayType{ComplexType}(undef, (L,L,2))

    for i in 1:L
        prethermalize(state, fft_temp, m², L^3)
        @show i
        flush(stdout)
        save_state("/home/jkott/perm/modelH_2D/thermalized/thermalized_L_$(L)_mass_$(m²)_id_$(seed).jld2", state, m², i)
    end
end

main()
