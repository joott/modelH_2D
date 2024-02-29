cd(@__DIR__)

using JLD2
using CodecZlib

include("src/modelH.jl")

@init_state
arrays = make_temp_arrays(state)

@time thermalize(state, arrays, m², 1000)
@time thermalize(state, arrays, m², 1000)
@time thermalize(state, arrays, m², 1000)
