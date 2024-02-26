cd(@__DIR__)

using JLD2
using CodecZlib

include("src/modelH.jl")

function main()
    @init_state
    # fft_temp = ArrayType{ComplexType}(undef, (L,L,L,3))
    arrays = make_temp_arrays(state)

    for i in 1:L^2
        thermalize(state, arrays, m², L^2)
        @show i
        save_state("/home/shared/modelH/thermalized/thermalized_L_$(L)_id_$(seed).jld2", state, m²)
    end
end

main()
