using ProgressMeter
using DelimitedFiles

N = 80000
N_bin = 500


# =========== LCG ============
x = 57378 
function LCG(a::Int64,c::Int64,m::Int64)::Float64
    global x
    x=(a*x+c)%m
    return x/m
end


rnd_lst=Vector{Float64}()
@showprogress  for i in 1:N
    append!(rnd_lst,LCG(1664525,1013904233,2^32))
end
writedlm("LCG_julia.dat", rnd_lst)

# =========== Metropolis-Hasting ============
e=MathConstants.e
eps = 4.5
pdf = x->0.204292*(exp(-(x+2)^2/pi)+0.6*exp(-(x-2)^2/e))
N_acpt = 0
rnd_lst =Vector{Float64}()
x = 0
@showprogress for i in 1:N
    x_cndt = x + (1-2*rand())*eps
    ds = pdf(x_cndt)/pdf(x)
    if pdf(x_cndt)/pdf(x) > rand()
        global N_acpt = N_acpt +1 
        global x = x_cndt
    end
    append!(rnd_lst,x)
end
println("accptnc: ", N_acpt/N)

writedlm("Mtrpls_Hstng_julia.dat", rnd_lst)
