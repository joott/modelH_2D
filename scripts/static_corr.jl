cd("/home/jkott/perm/modelH_2D/thermalized")
using DelimitedFiles
using Random
using Glob
using JLD2
using CodecZlib
using FFTW

const L = parse(Int, ARGS[1])

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

function collect_files()
    dfs = glob("thermalized_mass_-3.*_L_$(L)_*")
    println(length(dfs))
    collection = Dict{Float64, Vector{Array{Float64}}}()
    for filename in dfs
        file = load(filename)
        phi = file["u"][:,:,3]
        phik = fft(phi)
        phik2 = phik .* conj.(phik)
        G = ifft(phik2)[:,1]
        mass = file["m²"]

        if haskey(collection, mass)
            push!(collection[mass], real.(G))
        else
            collection[mass] = [real.(G)]
        end
    end

    return collection
end

function main()
    c = collect_files()

    println(keys(c))
    for mass in keys(c)
        mass_round = round(mass, digits=3)
        jldopen("static_H_$(L)_$(mass_round).jld2", "w") do savefile
            datasets = c[mass]
            N = length(datasets)
            println("for $(mass_round) we have $N")
            G = average(datasets)
            err = sqrt.(variance(datasets)/N)
            x = [0:L-1;] / L

            savefile["x"] = x
            savefile["G"] = G
            savefile["err"] = err
            savefile["N"] = N
            savefile["mass"] = mass_round
        end
    end
end

main()
