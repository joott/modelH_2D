using ParallelStencil
using ParallelStencil.FiniteDifferences2D

if cpu
    using FFTW
else
    using CUDA.CUFFT
end

include("helpers.jl")

const plans = (plan_fft(ArrayType{FloatType}(undef,(L,L,2)), (1,2)),
               plan_ifft!(ArrayType{ComplexType}(undef,(L,L,2)), (1,2)))

@parallel_indices (ix,iy) function poisson_scaling(πfft)
    k1 = sin(2pi * (ix-1) / L)
    k2 = sin(2pi * (iy-1) / L)
    k_factor = k1^2 + k2^2

    if k_factor > 1e-11
        k_factor = (k1 * πfft[ix,iy,1] + k2 * πfft[ix,iy,2]) / k_factor
        πfft[ix,iy,1] -= k1 * k_factor
        πfft[ix,iy,2] -= k2 * k_factor
    else
        πfft[ix,iy,1] = 0.0
        πfft[ix,iy,2] = 0.0
    end

    return
end

function project(π, temp)
    temp .= plans[1] * π
    @parallel (1:L, 1:L) poisson_scaling(temp)
    plans[2] * temp
    π .= real.(temp)

    !cpu && CUDA.reclaim()
end

##
@static if H0
  @static if NModC 
    @parallel function deterministic_elementary_step(
        π1, π2, ϕ,
        dπ1, dπ2, dϕ)

      return
    end
  else
    @parallel function deterministic_elementary_step(
        π1, π2, ϕ,
        dπ1, dπ2, dϕ)

      ### phi update
      # π_μ ∇_μ ϕ
      @all(dϕ) = -1.0/ρ * (@all(π1) * @d_xc(ϕ) + @all(π2) * @d_yc(ϕ))

      ### pi update
      # ∇_μ ϕ ∇²ϕ
      @all(dπ1) = -@d_xc(ϕ) * @d2_xy(ϕ)
      @all(dπ2) = -@d_yc(ϕ) * @d2_xy(ϕ)

      return
    end
  end 
else
  @static if NModC 
    @parallel function deterministic_elementary_step(
        π1, π2, ϕ,
        dπ1, dπ2, dϕ)

      # π_ν ∇_ν π_μ
      @all(dπ1) = -0.5/ρ * (@all(π1) * @d_xc(π1) + @all(π2) * @d_yc(π1))
      @all(dπ2) = -0.5/ρ * (@all(π1) * @d_xc(π2) + @all(π2) * @d_yc(π2))

      # ∇_ν π_μ π_ν
      @all(dπ1) = @all(dπ1) - 0.5/ρ * (@prd_d_xc(π1,π1) + @prd_d_yc(π1,π2))
      @all(dπ2) = @all(dπ2) - 0.5/ρ * (@prd_d_xc(π2,π1) + @prd_d_yc(π2,π2))

      return
    end
  else 
    @parallel function deterministic_elementary_step(
        π1, π2, ϕ,
        dπ1, dπ2, dϕ)

      ### phi update
      # π_μ ∇_μ ϕ
      @all(dϕ) = -1.0/ρ * (@all(π1) * @d_xc(ϕ) + @all(π2) * @d_yc(ϕ))

      ### pi update
      # ∇_μ ϕ ∇²ϕ
      @all(dπ1) = -@d_xc(ϕ) * @d2_xy(ϕ)
      @all(dπ2) = -@d_yc(ϕ) * @d2_xy(ϕ)

      # π_ν ∇_ν π_μ
      @all(dπ1) = @all(dπ1) - 0.5/ρ * (@all(π1) * @d_xc(π1) + @all(π2) * @d_yc(π1))
      @all(dπ2) = @all(dπ2) - 0.5/ρ * (@all(π1) * @d_xc(π2) + @all(π2) * @d_yc(π2))

      # ∇_ν π_μ π_ν
      @all(dπ1) = @all(dπ1) - 0.5/ρ * (@prd_d_xc(π1,π1) + @prd_d_yc(π1,π2))
      @all(dπ2) = @all(dπ2) - 0.5/ρ * (@prd_d_xc(π2,π1) + @prd_d_yc(π2,π2))

      return
    end 
  end

end
##
#
#


function deterministic(state, k1, k2, k3, rk_state, fft_temp)
    project(state.π, fft_temp)
    @parallel deterministic_elementary_step(view_tuple(state.u)..., view_tuple(k1)...)

    rk_state.u .= state.u .+ Δtdet*k1
    project(rk_state.π, fft_temp)
    @parallel deterministic_elementary_step(view_tuple(rk_state.u)..., view_tuple(k2)...)

    rk_state.u .= state.u .+ Δtdet*0.25*(k1 .+ k2)
    project(rk_state.π, fft_temp)
    @parallel deterministic_elementary_step(view_tuple(rk_state.u)..., view_tuple(k3)...)

    state.u .+= Δtdet*(0.5*k1 .+ 0.5*k2 .+ 2.0*k3)/3.0  
    project(state.π, fft_temp)
end

"""
  Elementary stochastic step with the transfer of the momentum density (μ-th component) from the cell x1 to x2 
"""
function pi_step(π, n, m, μ, (i,j))
    xy = ((2i + m)%L+1, j%L+1)
    x1 = (xy[n+1], xy[2-n])
    x2 = ((x1[1]-n)%L+1, (x1[2]-1+n)%L+1)

    norm = cos(2pi*rand())*sqrt(-2.0*log(rand()))
    q = Rate_pi * norm

    δH = (q * (π[x1..., μ] - π[x2..., μ]) + q^2)/ρ
    P = exp(-δH)
    r = rand()

    π[x1..., μ] += q * (r<P)
    π[x2..., μ] -= q * (r<P)
end

"""
  Computing the local change of energy in the cell x 
"""
function ΔH_phi(ϕ, m², x, q)
    ϕold = ϕ[x...]
    ϕt = ϕold + q
    Δϕ = ϕt - ϕold
    Δϕ² = ϕt^2 - ϕold^2

    ∑nn = (ϕ[NNp(x[1]), x[2]] + ϕ[x[1], NNp(x[2])]
         + ϕ[NNm(x[1]), x[2]] + ϕ[x[1], NNm(x[2])])

    return 2Δϕ² - Δϕ * ∑nn + 0.5m² * Δϕ² + 0.25λ * (ϕt^4 - ϕold^4)
end

function phi_step(ϕ, m², n, m, (i,j))
    xy = ((4i + 2j + m%2)%L+1, (j + m÷2)%L+1)
    x1 = (xy[n+1], xy[2-n])
    x2 = ((x1[1]-n)%L+1, (x1[2]-1+n)%L+1)

    norm = cos(2pi*rand())*sqrt(-2*log(rand()))
    q = Rate_phi * norm

    δH = ΔH_phi(ϕ, m², x1, q) + ΔH_phi(ϕ, m², x2, -q) + q^2
    P = exp(-δH)
    r = rand()

    ϕ[x1...] += q * (r<P)
    ϕ[x2...] -= q * (r<P)
end

##
@static if cpu

function pi_sweep(π, n, m)
    Threads.@threads for l in 0:L^2-1
        μ = l ÷ (L^2÷2) + 1
        i = (l ÷ L) % (L ÷ 2)
        j = l % L

        pi_step(π, n, m, μ, (i,j))
    end
end

function phi_sweep(ϕ, m², n, m)
    Threads.@threads for l in 0:L^2÷4-1
        i = l ÷ L
        j = l % L

        phi_step(ϕ, m², n, m, (i,j))
    end
end

else

function _pi_sweep(π, n, m)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1
    stride = gridDim().x * blockDim().x

    for l in index:stride:L^2-1
        μ = l ÷ (L^2÷2) + 1
        i = (l ÷ L) % (L ÷ 2)
        j = l % L

        pi_step(π, n, m, μ, (i,j))
    end
end

function _phi_sweep(ϕ, m², n, m)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1
    stride = gridDim().x * blockDim().x

    for l in index:stride:L^2÷4-1
        i = l ÷ L
        j = l % L

        phi_step(ϕ, m², n, m, (i,j))
    end
end

_pi_sweep_temp  = @cuda launch=false _pi_sweep(CuArray{FloatType}(undef,(L,L,2)), 0, 0)
_phi_sweep_temp = @cuda launch=false _phi_sweep(CuArray{FloatType}(undef,(L,L)), zero(FloatType), 0, 0)

const N_pi = L^2÷2
config = launch_configuration(_pi_sweep_temp.fun)
const threads_pi = min(N_pi, config.threads)
const blocks_pi = cld(N_pi, threads_pi)

const N_phi = L^2÷4
config = launch_configuration(_phi_sweep_temp.fun)
const threads_phi = min(N_phi, config.threads)
const blocks_phi = cld(N_phi, threads_phi)

pi_sweep  = (π, n, m) -> _pi_sweep_temp(π, n, m; threads=threads_pi, blocks=blocks_pi)
phi_sweep = (ϕ, m², n, m) -> _phi_sweep_temp(ϕ, m², n, m; threads=threads_phi, blocks=blocks_phi)

end
##

function dissipative(state, m²)
    # pi update
    for n in 0:1, m in 0:1
        pi_sweep(state.π, n, m)
    end

    # phi update
    for n in 0:1, m in 0:3
        phi_sweep(state.ϕ, m², n, m)
    end
end

"""
    sum_check(x)

Checks if any x, or any of its entries are infinite or NAN 
"""
function sum_check(x)
    s = sum(x)
    isnan(s) || !isfinite(s)
end

function thermalize(state, arrays, m², N)
    for _ in 1:N
        if sum_check(state.ϕ) 
            break
        end
        dissipative(state, m²)

        deterministic(state, arrays...)
    end
end

function prethermalize(state, fft_temp, m², N)
    for _ in 1:N
        if sum_check(state.ϕ) 
            break
        end
        dissipative(state, m²)
        project(state.π, fft_temp)
    end
end
