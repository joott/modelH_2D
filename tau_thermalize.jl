cd(@__DIR__)

using JLD2
using CodecZlib

include("src/modelH.jl")

function main()
    @init_state
    fft_temp = ArrayType{ComplexType}(undef, (L,L,L,3))

    prethermalize(state, fft_temp, m², L^3)
    tau_id = round(Int, 100 * parsed_args["mass"])

    for i in 1:L
        prethermalize(state, fft_temp, m², L^3)
        @show i
        flush(stdout)
        save_state("/home/jkott/perm/modelH/thermalized/tau/$(tau_id)/thermalized_L_$(L)_id_$(seed).jld2", state, m²)
    end
end

main()
