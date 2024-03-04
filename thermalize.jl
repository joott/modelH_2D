cd(@__DIR__)

using JLD2
using CodecZlib

include("src/modelH.jl")

function main()
    @init_state
    fft_temp = ArrayType{ComplexType}(undef, (L,L,L,3))

    for i in 1:L^2
        prethermalize(state, fft_temp, m², L^2)
        @show i
        save_state("/home/shared/modelH_2D/thermalized/thermalized_L_$(L)_id_$(seed).jld2", state, m²)
    end
end

main()
