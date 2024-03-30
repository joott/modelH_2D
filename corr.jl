cd(@__DIR__)

using DelimitedFiles
using JLD2
using CodecZlib

const N = 64
const mass = ARGS[1]

function autocor_loc_2(x, beg, max, n=2)
	C = zeros(Complex{Float64},max+1)
	N = zeros(Int64,max+1)
	Threads.@threads for tau in 0:max
		for i in beg:length(x)-max
			j = i + tau
			@inbounds C[tau+1] = C[tau+1] +  (x[i]*conj(x[j]))^n
			@inbounds N[tau+1] = N[tau+1] + 1
		end
	end
	(collect(0:max),  C ./ N)
end

df=readdlm("output_H_phi_L_$(N)_eta_0.1_mass_$(mass).dat",' ')
(t, C) = autocor_loc_2((df[:,2].+df[:,3].*1.0im), 1, N^3, 1)

jldsave("corr_phi_L_$(N)_mass_$(mass).dat", true; t, C)
