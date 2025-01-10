using ArgParse
using Distributions
using Random
using JLD2
using CodecZlib
using ParallelStencil

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--mass"
            help = "mass shift from m²c"
            arg_type = Float64
            default = 0.0
        "--dt"
            help = "size of time step"
            arg_type = Float64
            default = 0.04
        "--rho"
            help = "mass density parameter; defaults to 1.0"
            arg_type = Float64
            default = 1.0
        "--rng"
            help = "seed for random number generation"
            arg_type = Int
            default = 0
        "--init"
            help = "path of .jld2 file with initial state"
            arg_type = String
        "--fp64"
            help = "flag to use Float64 type rather than Float32"
            action = :store_true
        "--cpu"
            help = "parallelize on CPU rather than GPU"
            action = :store_true
        "--H0"
            help = "disable ππ term in deterministic step"
            action = :store_true
        "--NModC"
            help = "disable π-phi coupling in pi-deterministic step"
            action = :store_true
        "size"
            help = "side length of lattice"
            arg_type = Int
            required = true
        "viscosity"
            help = "the shear viscosity η in lattice units"
            arg_type = Float64
            required = true
    end

    return parse_args(s)
end

#=
 Parameters below are
 1. L is the number of lattice sites in each dimension; it accepts the second argument passed to julia   
 2. λ is the 4 field coupling
 3. Γ is the scalar field diffusion rate; in our calculations we set it to 1, assuming that the time is measured in the appropriate units 
 4. T is the temperature 
 5. m² = -2.28587 is the critical value of the mass parameter 
=#

parsed_args = parse_commandline()

const cpu = parsed_args["cpu"]
!cpu && using CUDA

const H0 = parsed_args["H0"]
const NModC = parsed_args["NModC"]
const FloatType = parsed_args["fp64"] ? Float64 : Float32
const ComplexType = complex(FloatType)
const ArrayType = cpu ? Array : CuArray
const SubArrayType = cpu ? SubArray : CuArray

const λ = FloatType(4.0)
const Γ = FloatType(1.0)
const T = FloatType(1.0)

const L = parsed_args["size"]
const η = FloatType(parsed_args["viscosity"])
const ρ = FloatType(parsed_args["rho"])
const m²c = FloatType(-3.824)
const m² = FloatType(m²c + parsed_args["mass"])
const Δt = FloatType(parsed_args["dt"]/Γ)

const Δtdet = Δt
const Rate_phi = FloatType(sqrt(2.0*Δt*Γ))
const Rate_pi  = FloatType(sqrt(2.0*Δt*η))
const ξ = Normal(FloatType(0.0), FloatType(1.0))

const seed = parsed_args["rng"]
if seed != 0
    Random.seed!(seed)
    !cpu && CUDA.seed!(seed)
end

struct State
    u::ArrayType
    π::SubArrayType
    ϕ::SubArrayType
    State(u) = new(u, @view(u[:,:,1:2]), @view(u[:,:,3]))
end

function hotstart(n, n_components)
	u = rand(ξ, n, n, n_components)

    for i in 1:n_components
        u[:,:,i] .-= shuffle(u[:,:,i])
    end

    State(ArrayType(u))
end

init_arg = parsed_args["init"]

##
if isnothing(init_arg)

macro init_state() esc(:( state = hotstart(L, 3) )) end

else

macro init_state()
    file = jldopen(init_arg, "r")
    state = State(ArrayType(file["u"]))
    return esc(:( state = $state ))
end

end
##

##
@static if cpu

@init_parallel_stencil(Threads, FloatType, 2);

else

@init_parallel_stencil(CUDA, FloatType, 2);

end
##
