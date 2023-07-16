using HDF5
using QuadGK
using GSL


dQ = parse(Float64,ARGS[1])
Q_end = parse(Float64,ARGS[2])

#? ==== running αs ====
αs_h5 =  h5open("./alphas.h5","r")
scale_ary = Array(αs_h5["/Q"])
αs_ary_LO = Array(αs_h5["/nloop_1/alphas"])
NQ = length(scale_ary)

αs_intrp = interp_alloc(gsl_interp_cspline, NQ)
acc_αs_intrp = interp_accel_alloc()
interp_init(αs_intrp, scale_ary, αs_ary_LO, NQ)

αs(Q::Float64)::Float64 = interp_eval(αs_intrp, scale_ary, αs_ary_LO, Q, acc_αs_intrp)


x_ary = [0.0112884, 0.0162378, 0.0233572, 0.0335982, 0.0483293, 0.0695193,0.1, 0.101, 0.148263, 0.195526, 0.242789, 0.290053, 0.337316,0.384579, 0.431842, 0.479105, 0.526368, 0.573632, 0.620895, 0.668158,0.715421, 0.762684, 0.809947, 0.857211, 0.904474, 0.951737, 0.999]

#? μ=sqrt(2)
#Q_strt = sqrt(2.0)
#pdf_ary = [1.69915, 1.91549, 2.11076, 2.34758, 2.54498, 2.68931, 2.70135, 2.6995, 2.50072, 2.20909, 1.88426, 1.56784, 1.27211, 1.02221, 0.80981, 0.636259, 0.48584, 0.36576, 0.276329, 0.211289, 0.149932, 0.106993, 0.0792276, 0.0505565, 0.0267606, 0.0029888, 0.00622659]

#? μ=2
Q_strt = 2.0
pdf_ary = [1.16823, 1.23051, 1.32372, 1.44673, 1.59438, 1.75377, 1.89938,1.90281, 1.9919, 1.97738, 1.8986, 1.77894, 1.63264, 1.46887, 1.29438,1.11501, 0.936497, 0.764578, 0.604703, 0.461515, 0.33831, 0.23663,0.156137, 0.0948638, 0.0498701, 0.0183754, 0.000166569]

append!(x_ary,1.0)
append!(pdf_ary,0.0)

Nx = length(x_ary)

pdf_intrp = interp_alloc(gsl_interp_cspline, Nx)
acc_pdf_intrp = interp_accel_alloc()


function dpdfdQ(intrp,Q::Float64)::Vector{Float64}
    pdf(x::Float64) = interp_eval(intrp, x_ary, pdf_ary, x, acc_pdf_intrp)
    #intgrnd = z::Float64-> (z*(pdf.(x_ary/z) - z*pdf.(x_ary))/(1-z))
    intgrl_ary = zeros(Nx)
    for i in 1:Nx-1
        intgrnd = z::Float64-> (pdf(x_ary[i]/z) - pdf(x_ary[i]))/(1-z)
        intgrl_ary[i] = quadgk(intgrnd, x_ary[i], 0.999, rtol=1e-6)[1]
    end
    #res =  αs(Q)/π/Q*( (2.0 .+ 8.0/3*log.(1 .- x_ary)) .* pdf.(x_ary) + 8.0*intgrl_ary/3 )
    res = @. αs(Q)/π/Q*( (2.0 + 8.0/3*log(1 - x_ary)) * pdf(x_ary) + 8.0*intgrl_ary/3 )
    res[Nx] = 0.0
    return res
end


Q_ary = [Q_strt]
pdf_evlvd_ary = [pdf_ary]
Q = Q_strt 
while Q_strt<=Q<=Q_end    

    #? ==== K1 = dQ*RHS(Q,pdf) =====
    interp_init(pdf_intrp, x_ary, pdf_ary, Nx)
    K1_ary = dQ * dpdfdQ(pdf_intrp, Q)
    #println(" === K1 done ! ===")
    
    #? ==== K2 = dQ*RHS(Q+dQ/2, pdf + K1/2) =====
    K2_pdf_ary = pdf_ary + K1_ary/2
    interp_init(pdf_intrp, x_ary, K2_pdf_ary, Nx)
    K2_ary = dQ * dpdfdQ(pdf_intrp, Q+dQ/2)
    #println(" === K2 done ! ===")

    #? ==== K3 = dQ*RHS(Q+dQ/2, pdf + K2/2) =====
    K3_pdf_ary = pdf_ary + K2_ary/2
    interp_init(pdf_intrp, x_ary, K3_pdf_ary, Nx)
    K3_ary = dQ * dpdfdQ(pdf_intrp,Q+dQ/2)
    #println(" === K3 done ! ===")

    #? ==== K4 = dQ*RHS(Q+dQ, pdf + K3) =====
    K4_pdf_ary = pdf_ary + K3_ary
    interp_init(pdf_intrp, x_ary, K4_pdf_ary, Nx)
    K4_ary = dQ * dpdfdQ(pdf_intrp,Q+dQ)
    #println(" === K4 done ! ===")

    global Q = Q + dQ
    global pdf_ary = pdf_ary + K1_ary/6 + K2_ary/3 + K3_ary/3 + K4_ary/6
    global pdf_ary[Nx] = 0.0

    append!(Q_ary,Q)
    append!(pdf_evlvd_ary,[pdf_ary])

    println(" === evolve to $Q done ! ===")
    interp_init(pdf_intrp, x_ary, pdf_ary, Nx)

end

evltn_file=h5open("./trans_evltn.h5","w")
write(evltn_file,"/non_singlet/trnsvrsty",hcat(pdf_evlvd_ary...))
write(evltn_file,"/Q",hcat(Q_ary...))
write(evltn_file,"/x",hcat(x_ary...))
close(evltn_file)