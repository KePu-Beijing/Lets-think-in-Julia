using ProgressMeter
using DataFrames, Gadfly, Cairo, Fontconfig

N_t = parse(Int64, ENV["Nt"]) #20
a = parse(Float64, ENV["a"]) #0.5
N_path = parse(Int64, ENV["Npath"]) #20000
ϵ = parse(Float64, ENV["a"]) #1.0
N_corr = parse(Int64, ENV["Ncorr"]) #20

N_thrml = 4*N_corr;

function S(x::Vector{Float64}, i::Int64)::Float64
    return a*x[i]^2/2 +x[i]*(x[i]-x[mod(i+1-1, N_t)+1]-x[mod(i-1-1, N_t)+1])/a
end

N_acpt = 0
function update!(x::Vector{Float64})
    for i in 1:N_t
        x_i= x[i]
        S_i = S(x,i)
        x[i] = x[i] + (1-2*rand())*ϵ
        global N_acpt = N_acpt + 1
        dS = S(x, i) - S_i
        if dS>0 && exp(-dS) < rand()
            x[i] = x_i
            global N_acpt = N_acpt - 1
        end
    end
end

x = zeros(N_t)
N_acpt = 0
println(">>> Thermalizing ......")
@showprogress for i in 1:N_thrml
    update!(x)
end
println("    Acceptance: ", N_acpt/(N_thrml*N_t))

N_acpt = 0
x_lst = []
println(">>> Gnerating paths ......")
@showprogress for i in 1:N_path
        for j in 1:N_corr
            update!(x)
        end
        append!(x_lst,[Vector(x)])        
    end


println(">>> Calculating 2-point function ......")



function C2_func_p(x::Vector{Float64},n::Int64)::Float64
    return x[1]*x[(1-1+n)%N_t+1]
end

function C2_func(x::Vector{Float64},n::Int64)::Float64
    return sum(1:N_t .|> (i->x[i]*x[(i-1+n)%N_t+1]))/N_t
end

C2_lst = [0:N_t-1 .|> (n->C2_func(x_lst[i],n)) for i in 1:N_path]
C2_avg = sum(C2_lst)/N_path
C2_err = sqrt.(sum([(C2_lst[i] - C2_avg).^2/N_path for i in 1:N_path])/N_path)


N_t_dE = N_t ÷ 3 +1

dE_avg = 1/a * log.(C2_avg[1:N_t_dE-1] ./ C2_avg[2:N_t_dE])
dE_err = 1/a * sqrt.((C2_err[1:N_t_dE-1] ./ C2_avg[1:N_t_dE-1]).^2+(C2_err[2:N_t_dE] ./ C2_avg[2:N_t_dE]).^2);

C2_df = DataFrame(x=0:length(C2_avg)-1, y=C2_avg,mins=C2_avg-C2_err, maxs=C2_avg+C2_err);
dE_df = DataFrame(x=0:N_t_dE-1-1, y=dE_avg,mins=dE_avg-dE_err, maxs=dE_avg+dE_err);

println(">>> Plotting ...... ")

C2_plot=plot(C2_df, x=0:length(C2_avg)-1, y=:y, ymin=:mins, ymax=:maxs,Geom.point, Geom.errorbar,size=[3pt], Scale.y_log10,Theme(background_color="white",default_color="midnightblue"));
dE_plot=plot(dE_df, x=0:N_t_dE-1-1, y=:y, ymin=:mins, ymax=:maxs,yintercept=[1], Geom.hline(color="green",style=:dot),Geom.point, Geom.errorbar,size=[3pt],
             Coord.cartesian(xmin=-0.2 ,ymin=0, ymax=2),Theme(background_color="white",default_color="midnightblue"));
set_default_plot_size(1200px, 400px)
plots=hstack(C2_plot,dE_plot)

draw(PNG("hrmnc_oscltr_julia.png", 1200px, 400px, dpi=400), plots)
