cd(@__DIR__)

using JLD2
using CodecZlib
using Printf

include("src/modelH.jl")

function op(phi)
    phik = cpu ? fft(phi) : Array(CUFFT.fft(phi))
    phik[:,1,1]
end

function main()
    @init_state
    arrays = make_temp_arrays(state)

    thermalize(state, arrays, m², L^4)

    eta = round(η,digits=3)
    model = H0 ? (NModC ? "diff" : "H0") : (NModC ? "self" : "H")

    skip = L
    maxt = 4*L^4

    open("/home/jkott/perm/modelH_2D/dynamics/output_$(model)_phi_L_$(L)_eta_$(eta)_mass_$(m²).dat", "w") do iophi
        for idx in 0:skip:maxt
            ϕk = op(state.ϕ)
            Printf.@printf(iophi, "%i", idx)
            for kx in 2:4
                Printf.@printf(iophi, " %.15f %.15f", real(ϕk[kx]), imag(ϕk[kx]))
            end

            Printf.@printf(iophi, "\n")

            if idx%L^4 == 0
                save_state("/home/jkott/perm/modelH_2D/thermalized/thermalized_L_$(L)_mass_$(m²)_id_$(seed).jld2", state, m², idx)
            end

            thermalize(state, arrays, m², skip)
        end
    end
end

main()
