mom_ary = [0.0, 0.545, ... ,17.44]

wrkng_dir = pwd()

for P in mom_ary
  P_id = replace(string(P), "." => "d")

  @eval Γ_plus(n::Int64, k::Float64) = ccall(($(string("Vrtx_", P_id, "_plus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)
  @eval Γ_minus(n::Int64, k::Float64) = ccall(($(string("Vrtx_", P_id, "_minus")), "$wrkng_dir/libGPD_3meson_vrtx_$(prt).so"), Cdouble, (Cint, Cdouble), n, k)


  
  ...

end