using HDF5
using QuadGK
using GSL
using ProgressMeter

#!!!indx start from 0, as it is pass to C function in .so

indx = parse(Int, ARGS[1])
nA = parse(Int, ARGS[2])
prt = ARGS[3]

#P_array = [2.18, 4.36, 8.72, 13.08, 17.44]
#ξ_array = [0.25, 0.5, 0.75]

N_sum = 47

mom_ary = [1.12875, 2.2575, 3.38625, 4.515, 5.64375, 6.7725, 7.90125, 9.03, 11.2875, 13.545, 15.8025, 18.06, 22.575, 27.09, 31.605, 36.12, 45.15, 54.18, 63.21]

#=
Γ_vrtx_ids = ["P2d18_xi0d25", "P2d18_xi0d5",   "P2d18_xi0d75",  "P4d36_xi0d25", "P4d36_xi0d5",
              "P4d36_xi0d75", "P8d72_xi0d25",  "P8d72_xi0d5",   "P8d72_xi0d75", "P13d08_xi0d25",
              "P13d08_xi0d5", "P13d08_xi0d75", "P17d44_xi0d25", "P17d44_xi0d5", "P17d44_xi0d75"]
=#

include("./Gamma_vrtx_ids.txt")

findnearest = mom->mom_ary[findmin(abs.(mom_ary .- mom))[2]] #!!! get rid of rounding like 7.630000000000001~7.64  1.0899999999999999~1.09

function Pξ_parser(Γid)
    mtchd = match(r"^P(.*)_xi(.*)", Γid).captures
    return @. parse(Float64, replace(mtchd, "d" => "."))
end



#!!! BEAWARE P_indx is passed to GammaVrtx called from libGPD_3meson_vrtx.so, indx starts from 0, indicating the index of momentum, xi combinations
#!!! see def in libGPD_3meson_vrtx.cpp, double (*ptr_GammaVrtx_k12_intgrnd_P[N_moms])(int nA, int nC, double k1, double k2) = ...
#!!! here in julia, index starts from 1

P_val, ξ = Pξ_parser(Γ_vrtx_ids[indx+1])
Δ_val = 2 * ξ * P_val

P = findnearest(P_val)
Δ = findnearest(Δ_val)


wrkng_dir = pwd()


μ_ary = h5read("./h5_data/tHooft_WF_coeff_bsw_cc.h5", "/cc/N_384/mass")


P0 = nC -> sqrt(μ_ary[nC+1]^2 + P^2)
Δ0 = nC -> sqrt(μ_ary[nC+1]^2 + Δ^2)

P_id = replace(string(P), "." => "d")
Δ_id = replace(string(Δ), "." => "d")

ϵ = 1e-6

#@eval Γ_vrtx_k12(nA::Int64,nB::Int64,nC::Int64,k1::Float64,k2::Float64)=ccall(($(string("GammaVrtx_k12_intgrnd_P",P_id)),"$wrkng_dir/libGPD_3meson_vrtx.so"),Cdouble,(Cint,Cint,Cint,Cdouble,Cdouble),nA,nB,nC,k1,k2)        



@eval ψ_plus(n::Int64, k::Float64) = ccall(($(string("BGWFP", P_id, "plus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)
@eval ψ_minus(n::Int64, k::Float64) = ccall(($(string("BGWFP", P_id, "minus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)

@eval ψ_plus_Δ(n::Int64, k::Float64) = ccall(($(string("BGWFP", Δ_id, "plus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)
@eval ψ_minus_Δ(n::Int64, k::Float64) = ccall(($(string("BGWFP", Δ_id, "minus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)


#!!!indx start from 0, as it is pass to C function in .so
@eval GammaVrtx(nA::Int64, nC::Int64, P_indx::Int64) = ccall($(("GammaVrtx_$prt"), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cint, Cint), nA, nC, P_indx)

@eval θ(k::Float64) = ccall(("theta", "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cdouble,), k)


Λ = 40

#(p-P,-P) -> (p,P) => (p+P -P, -P) = (p,-P)  -> (p+P,P) = (-1)^n (-p,P)
#(p,-P) -> (-1)^n (-p,P)


#!! without sqrt(2q_nc^2) in Zhewen Mo's note
function Ψ(nC::Int64, x::Float64, part::String)::Float64
    if part == "Ket"
        return (-1)^nC / (4π) * cos(θ(P * (x + ξ)) / 2 + θ(P * (x - ξ)) / 2) * (ψ_plus_Δ(nC, -(x - ξ) * P) - ψ_minus_Δ(nC, -(x - ξ) * P))
    elseif part == "Bra"
        return       1 / (4π) * cos(θ(P * (x + ξ)) / 2 + θ(P * (x - ξ)) / 2) * (ψ_plus_Δ(nC,  (x + ξ) * P) - ψ_minus_Δ(nC,  (x + ξ) * P))
    end
end


x_ary = Array(-1.6:0.02:1.6)

Nx = length(x_ary)


Γ_res = zeros(N_sum + 1)
Ψ_res = zeros(N_sum + 1, Nx)

m = μ_ary[nA+1]

p = Progress(N_sum + 1)
Threads.@threads for nC in 0:N_sum
    #!!! uncomment this to calc Gamma_vrtx
    #if isodd(nC)
        Γ_res[nC+1] = GammaVrtx(nA, nC, indx)
    #else
    #    Γ_res[nC+1] = 0.0
    #end
    f_tmp = x -> Ψ(nC, x, prt)
    Ψ_res[nC+1, :] = f_tmp.(x_ary)
    next!(p)
end


#!!! uncomment this to write Gamma_vrtx

open("./nA=$(nA)_txt_data/off_$(Γ_vrtx_ids[indx+1])_GammaVrtx.txt", "a") do io
    println(io, "off_$(Γ_vrtx_ids[indx+1])_GammaVrtx_$(prt) = ", Γ_res)
end


open("./nA=$(nA)_txt_data/off_$(Γ_vrtx_ids[indx+1])_Psi.txt", "a") do io
    println(io, "off_$(Γ_vrtx_ids[indx+1])_Psi_$(prt) = ", Ψ_res)
end
