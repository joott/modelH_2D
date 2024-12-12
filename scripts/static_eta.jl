cd("/home/jkott/perm/modelH_2D/thermalized")
using DelimitedFiles
using Random
using Glob
using JLD2
using CodecZlib
using FFTW

const L = parse(Int, ARGS[1])
const eta = ARGS[2]

function variance(x)
    S=zeros(length(x[1]))
    mean = average(x)
    for i in eachindex(x)
      S .= S .+ (x[i] .- mean) .^2
    end
    S./length(x)
end

function average(x)
    sum(x)/length(x)
end

function bootstrap(fs, M)
    len = length(fs)
    bs_tot = Vector{Vector{Float64}}(undef, M)
    Threads.@threads for i in 1:M
        bsD = Vector{Vector{Float64}}(undef, len)
        for j in 1:len
            bsD[j] = fs[rand(1:len)]
        end
        bs_tot[i] = average(bsD)
    end
    mean = average(bs_tot)
    var = variance(bs_tot)
    (mean, sqrt.(var))
end

function collect_data()
    dfs = glob("Cum_modelH0_$(eta)_-0.035/thermalized_L_$(L)_*")
    N = length(dfs)
    datasets = Vector{Vector{Float64}}(undef, N)

    for i in 1:length(dfs)
        phi = load(dfs[i])["u"][:,:,3]
        phik = fft(phi)
        phik2 = phik .* conj.(phik)
        G = ifft(phik2)[:,1,1]
        datasets[i] = real.(G)
    end

    G = average(datasets)
    err = sqrt.(variance(datasets)/N)

    x = [0:L-1;] / L
    (x, G, err, N)
end

function main()
    (x, G, err, N) = collect_data()

    jldopen("static_H_$(L)_$(eta).jld2", "w") do savefile
        savefile["x"] = x
        savefile["G"] = G
        savefile["err"] = err
        savefile["N"] = N
    end
end

main()
