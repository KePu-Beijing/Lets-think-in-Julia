using HDF5
using QuadGK
using GSL
using DelimitedFiles
#using LsqFit
using IMinuit
using Distributed
using SharedArrays



addprocs(30)


@everywhere using Distributed, HDF5, QuadGK, DelimitedFiles, IMinuit, GSL, SharedArrays


@everywhere begin

    μ = parse(Float64, ENV["mu"])
    alttc = parse(Float64, ENV["alttc"])
    Pzlttc = parse(Int64, ENV["Pzlttc"])
    Smr = ENV["smr"]
    ΛQCD = parse(Float64, ENV["LmbdQCD"])
    ΛMS = parse(Float64, ENV["LmbdMS"])
    NL = parse(Int64, ENV["NL"])
    zs = parse(Float64, ENV["zs"])
    PnDatfile = ENV["PnDatfile"]
    P0Datfile = ENV["P0Datfile"]
    Outfile = ENV["Outfile"]
    k = parse(Float64, ENV["k"])
    dPDF = parse(Float64, ENV["dPDF"])
    δsys = parse(Float64, ENV["errsys"])

    λ_extrplt_strt = parse(Float64, ENV["lmbd_strt"])
    λ_extrplt_end = parse(Float64, ENV["lmbd_end"])
    λ_lng = parse(Float64, ENV["lmbd_lng"])

    yMax = parse(Float64, ENV["yMax"])
    xMax = parse(Float64, ENV["xMax"])
    xMin = parse(Float64, ENV["xMin"])


    Ns = parse(Int64, ENV["Ns"])
    Nx = parse(Int64, ENV["Nx"])
    UV = parse(Float64, ENV["QuadGKUV"])
    ϵ = parse(Float64, ENV["eps"])

    btstrp_id_strt = parse(Int64,ENV["btstrp_id_strt"])
    btstrp_id_end = parse(Int64,ENV["btstrp_id_end"])

    Nbtstrp = btstrp_id_end - btstrp_id_strt + 1



    γE = 0.5772156649015329
    expγE = exp(2 * γE)

    GeVfm = 0.1973269631
    Nf = 3
    b0 = 11 - 2.0 / 3 * Nf
    CF = 4.0 / 3
    αs = 2π / (b0 * log(μ / ΛMS))

    αs1 = 2π / (9 * log(μ / 0.24451721864451428))

    #m0all = 0.584466

    m0all = 0.9999085

    #function f1PDF(z::Folat64)::Float64

    f1PDFallx = [0.0, 0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.6, 0.6599999999999999, 0.72, 0.78, 0.84, 0.8999999999999999, 0.96, 1.02, 1.08, 1.14, 1.2, 1.26, 1.3199999999999998, 1.38, 1.44, 1.5, 1.56, 1.6199999999999999, 1.68, 1.74, 1.8, 1.86, 1.92, 1.98, 2.04, 2.1, 2.16, 2.22, 2.28, 2.34, 2.4, 2.46, 2.52, 2.58, 2.64, 2.7, 2.76, 2.82, 2.88, 2.94, 3.]
    f1PDFally = [0.358022, -0.195807, -0.433636, -0.488038, -0.508858, -0.512462, -0.512422, -0.49898, -0.468816, -0.438756, -0.414745, -0.382962, -0.351789, -0.285242, -0.217616, -0.188122, -0.117031, -0.0737374, -0.128576, -0.254514, -0.419045, -0.768752, -1.30363, -2.02369, -2.92893, -4.01934, -5.29492, -6.75568, -8.40162, -10.2327, -12.249, -14.4505, -16.8371, -19.4089, -22.1659, -25.1081, -28.2354, -31.548, -35.0457, -38.7285, -42.5966, -46.6498, -50.8882, -55.3118, -59.9205, -64.7144, -69.6935, -74.8578, -80.2073, -85.7419, -91.4617]
    #f1PDFally = [0.43025, -0.123819, -0.361709, -0.416307, -0.437171, -0.440704,-0.440611, -0.427088, -0.396668, -0.366311, -0.342009, -0.310077, -0.279048, -0.212751, -0.145287, -0.115959, -0.0452188, -0.00255544, -0.0580563, -0.18446, -0.349521, -0.699726, -1.23508, -1.95558, -2.86122, -3.95201, -5.22794, -6.68902, -8.33525, -10.1666, -12.1831, -14.3848, -16.7716, -19.3436, -22.1007, -25.0429, -28.1703, -31.4828, -34.9805, -38.6633, -42.5313, -46.5844, -50.8227, -55.2461, -59.8547, -64.6484, -69.6272, -74.7912, -80.1404, -85.6746, -91.3941]
    Nf1PDFall = length(f1PDFallx)
    f1PDFall_intrp = interp_alloc(gsl_interp_steffen, Nf1PDFall)
    acc_f1PDFall_intrp = interp_accel_alloc()
    interp_init(f1PDFall_intrp, f1PDFallx, f1PDFally, Nf1PDFall)
    f1PDFM(z::Float64)::Float64 = interp_eval(f1PDFall_intrp, f1PDFallx, f1PDFally, z, acc_f1PDFall_intrp)

    #println("== f1PDFM ==>",f1PDFM.([0.4,1.2]))

    function ZRPDF(a::Float64, z::Float64, Λ::Float64, k::Float64, m0::Float64, d::Float64)::Float64
        aG = a / GeVfm
        return exp(k * z / (aG * log(aG * Λ)) + m0 * z + f1PDFM(z) * aG^2 + (3 * CF) / b0 * log(log(1 / (aG * Λ)) / log(μ / Λ)) + log(1 + d / log(aG * Λ)))
    end


    Pz = Pzlttc * (2π) / (NL * alttc) * GeVfm
    Pzfile = h5open(PnDatfile, "r")
    
    
    #? P0 matrix elements
    P0file = h5open(P0Datfile, "r")
    P0C0Re = Array(P0file["/1/c0/Re"])
    P0C0Re = P0C0Re ./ P0C0Re[1]


    Nlnks = length(P0C0Re)
    lnks_Max = Nlnks - 1

    z_ary = Array(alttc .* (0:lnks_Max))

    #? λ = Pz ⋅ z, z in fm, Pz in fm^-1
    λ_ary = Pz .* z_ary / GeVfm


    λs = zs * Pz / GeVfm

    #println("=== λs ===>",[zs,Pz,GeVfm,λs])

    ZRPDF_ary = ZRPDF.(alttc, z_ary, ΛQCD, k, m0all, dPDF)

    GRmdl(λ::Float64, p::Vector{Float64}) = exp(-(λ / p[5])) * (p[1] * λ^(-p[3]) * cos((p[3] * π) / 2) + p[2] * λ^(-p[4]) * cos((p[4] * π) / 2 - λ))
    GImdl(λ::Float64, p::Vector{Float64}) = exp(-(λ / p[5])) * (-p[1] * λ^(-p[3]) * sin((p[3] * π) / 2) + p[2] * λ^(-p[4]) * sin((p[4] * π) / 2 - λ))


    λ_IR_ary = LinRange(λ_extrplt_strt, λ_extrplt_end, 41)

    Nλ_IR = length(λ_IR_ary)

    y_ary = Array(LinRange(-yMax, yMax, Ns))
    x_ary = Array(LinRange(xMin, xMax, Nx))
    filter!(x -> x != 0, x_ary)

    #? matching factor
    function Zh(x::Float64, Pz::Float64)
        if 0 < x < 1
            return αs1 * CF / (2π) * (((2x) / (1 - x)) * (log(4 * Pz^2 / μ^2) + log(x * (1 - x))) - (2x) / (1 - x) + 2 / (1 - x)) + (αs1 * CF / π) * (-1 / abs(1 - x) + 2 * (sf_Si((1 - x) * abs((zs / GeVfm) * Pz)) / (π * (1 - x))))
        elseif x > 1
            return αs1 * CF / (2π) * (((2x) / (1 - x)) * log(x / (x - 1)) - 2 / (1 - x)) + (αs1 * CF / π) * (-1 / abs(1 - x) + 2 * (sf_Si((1 - x) * abs((zs / GeVfm) * Pz)) / (π * (1 - x))))
        else
            return αs1 * CF / (2π) * (-((2x) / (1 - x)) * log(x / (x - 1)) + 2 / (1 - x)) + (αs1 * CF / π) * (-1 / abs(1 - x) + 2 * (sf_Si((1 - x) * abs((zs / GeVfm) * Pz)) / (π * (1 - x))))
        end
    end


    #? arry to store C0 read from all bootstrap samples
    c0Re_ary = [Array(Pzfile["/$(i)/c0/Re"]) for i in btstrp_id_strt:btstrp_id_end]
    c0Im_ary = [Array(Pzfile["/$(i)/c0/Im"]) for i in btstrp_id_strt:btstrp_id_end]

    
    #LCPDF_ary_shrd = zeros(Nbtstrp,length(x_ary)) 

end

#? shared arry to store bootstrap LCPDF after renomalize, FT, match quasi-PDF ...
LCPDF_ary_shrd=SharedMatrix{Float64}(Nbtstrp,length(x_ary))

@everywhere function LCPDF_proc(pid::Int64,LCPDF_ary_shrd)

    qPDF_ary = zeros(Ns)
    LCPDF_ary = zeros(length(x_ary))

    C0Re_slf_intrp = interp_alloc(gsl_interp_steffen, Nlnks)
    C0Im_slf_intrp = interp_alloc(gsl_interp_steffen, Nlnks)
    C0Re_hybrd_intrp = interp_alloc(gsl_interp_steffen, Nlnks)
    C0Im_hybrd_intrp = interp_alloc(gsl_interp_steffen, Nlnks)

    acc_C0Re_slf_intrp = interp_accel_alloc()
    acc_C0Im_slf_intrp = interp_accel_alloc()
    acc_C0Re_hybrd_intrp = interp_accel_alloc()
    acc_C0Im_hybrd_intrp = interp_accel_alloc()

    qPDF_intrp = interp_alloc(gsl_interp_steffen, Ns)
    acc_qPDF_intrp = interp_accel_alloc()

    #? bare matrix element C_0(z)
    #c0Re = Array(Pzfile["/$(i)/c0/Re"])
    #c0Im = Array(Pzfile["/$(i)/c0/Im"])

    c0Re = c0Re_ary[pid-btstrp_id_strt+1]
    c0Im = c0Im_ary[pid-btstrp_id_strt+1]

    #println(c0Re)

    c0Re0 = c0Re[1]
    c0Re = c0Re ./ c0Re0
    c0Im = c0Im ./ c0Re0

    C0Re_slf = c0Re ./ P0C0Re
    C0Im_slf = -1 .* c0Im ./ P0C0Re

    C0Re_hybrd = c0Re ./ ZRPDF_ary
    C0Im_hybrd = -1 .* c0Im ./ ZRPDF_ary

    interp_init(C0Re_slf_intrp, λ_ary, C0Re_slf, Nlnks)
    interp_init(C0Im_slf_intrp, λ_ary, C0Im_slf, Nlnks)

    interp_init(C0Re_hybrd_intrp, λ_ary, C0Re_hybrd, Nlnks)
    interp_init(C0Im_hybrd_intrp, λ_ary, C0Im_hybrd, Nlnks)

    #? renomalized matrix elements 

    # Interpolate
    C0Re_slf_F(λ::Float64)::Float64 = interp_eval(C0Re_slf_intrp, λ_ary, C0Re_slf, λ, acc_C0Re_slf_intrp)
    C0Im_slf_F(λ::Float64)::Float64 = interp_eval(C0Im_slf_intrp, λ_ary, C0Im_slf, λ, acc_C0Im_slf_intrp)

    C0Re_hybrd_F(λ::Float64)::Float64 = interp_eval(C0Re_hybrd_intrp, λ_ary, C0Re_hybrd, λ, acc_C0Re_hybrd_intrp)
    C0Im_hybrd_F(λ::Float64)::Float64 = interp_eval(C0Im_hybrd_intrp, λ_ary, C0Im_hybrd, λ, acc_C0Im_hybrd_intrp)


    ReC0_r = real((C0Re_slf_F(λs) + im * C0Im_slf_F(λs)) / (C0Re_hybrd_F(λs) + im * C0Im_hybrd_F(λs)))

    function C0Re_F(λ::Float64)::Float64
        if λ > λs
            return C0Re_hybrd_F(λ) * ReC0_r
        else
            return C0Re_slf_F(λ)
        end
    end

    function C0Im_F(λ::Float64)::Float64
        if λ > λs
            return C0Im_hybrd_F(λ) * ReC0_r
        else
            return C0Im_slf_F(λ)
        end
    end


    #? extrapolation 

    #C0Re_IR_ary_cmbnd = hcat(ones(Nλ_IR), λ_IR_ary, Array(C0Re_F.(λ_IR_ary)))
    #C0Im_IR_ary_cmbnd = hcat(2 * ones(Nλ_IR), λ_IR_ary, Array(C0Im_F.(λ_IR_ary)))

    #C0_IR_ary_cmbnd = vcat(C0Re_IR_ary_cmbnd, C0Im_IR_ary_cmbnd)


    VGRmdl(p::Vector{Float64}) = @. exp(-(λ_IR_ary / p[5])) * ( p[1] * λ_IR_ary^(-p[3]) * cos((p[3] * π) / 2) + p[2] * λ_IR_ary^(-p[4]) * cos((p[4] * π) / 2 - λ_IR_ary))
    VGImdl(p::Vector{Float64}) = @. exp(-(λ_IR_ary / p[5])) * (-p[1] * λ_IR_ary^(-p[3]) * sin((p[3] * π) / 2) + p[2] * λ_IR_ary^(-p[4]) * sin((p[4] * π) / 2 - λ_IR_ary))


    fcn(p) = sum((VGRmdl(p) - C0Re_F.(λ_IR_ary)).^2) + sum((VGImdl(p) - C0Im_F.(λ_IR_ary)).^2)
    pnames = ["p1", "p2", "p3", "p4", "p5"]
    #minuit = Minuit(fcn, [1,1,1,1,1]; name=pnames, error=1*ones(5),limit_p5=(0,10))
    minuit = Minuit(fcn, [1,1,1,2,5]; name=pnames, error=0.001*ones(5), limit_p3=(0,10),limit_p4=(0,10),limit_p5=(0,10))

    migrad(minuit)

    param_ary =  [minuit.values[i] for i in 1:5]
    
    param_ary #[CLSfit.args[1].args[2].args[i].args[2] for i in 1:5]

    function C0Re_extrpltd(λ::Float64)
        if λ > λ_lng
            return GRmdl(λ, param_ary)
        else
            return C0Re_F(λ)
        end
    end

    function C0Im_extrpltd(λ::Float64)
        if λ > λ_lng
            return GImdl(λ, param_ary)
        else
            return C0Im_F(λ)
        end
    end

    for i in 1:Ns
        y = y_ary[i]
        #global 
        intgrnd = λ::Float64 -> cos(λ * y) * C0Re_extrpltd(λ) - sin(λ * y) * C0Im_extrpltd(λ)

        qPDF_ary[i] = quadgk(intgrnd, 0, UV, rtol = 1e-6, atol = 1e-6)[1] / π

        #println(y_ary[i],"\t",qPDF_ary[i])
    end

    #? matching 
    interp_init(qPDF_intrp, y_ary, qPDF_ary, Ns)

    qPDF(y::Float64)::Float64 = interp_eval(qPDF_intrp, y_ary, qPDF_ary, y, acc_qPDF_intrp)

    for i in 1:length(x_ary)
        x = x_ary[i]
        #global 
        LCPDF_krnl(y::Float64)::Float64 = Zh(x / y, y * Pz) * qPDF(y) / abs(y) - Zh(y / x, x * Pz) * qPDF(x) / abs(x)

        x_brkpnts = [x - ϵ, -ϵ, ϵ, x + ϵ]
        sort!(x_brkpnts)

        if x >= 0
            LCPDF_ary[i] = qPDF(x) - (quadgk(LCPDF_krnl, -yMax, -ϵ, rtol = 1e-5, atol = 1e-5)[1] + quadgk(LCPDF_krnl, ϵ, x - ϵ, rtol = 1e-5, atol = 1e-5)[1] + quadgk(LCPDF_krnl, x + ϵ, yMax, rtol = 1e-5, atol = 1e-5)[1])
        else
            LCPDF_ary[i] = qPDF(x) - (quadgk(LCPDF_krnl, -yMax, x - ϵ, rtol = 1e-4, atol = 1e-4)[1] + quadgk(LCPDF_krnl, x + ϵ, -ϵ, rtol = 1e-4, atol = 1e-4)[1] + quadgk(LCPDF_krnl, ϵ, yMax, rtol = 1e-4, atol = 1e-4)[1])
        end

        # println(x_ary[i],"\t",LCPDF_ary[i])
    end

    #println("##### ",pid-btstrp_id_strt, " #####")

    LCPDF_ary_shrd[pid-btstrp_id_strt+1,:] = LCPDF_ary

    #println(" ===> $pid done !===>")


end

@time @sync @distributed for pid in btstrp_id_strt:btstrp_id_end
    LCPDF_proc(pid,LCPDF_ary_shrd)
end

#h5_file = Outfile

LCPDF_file=h5open(Outfile, isfile(Outfile) ? "r+" : "w")



for pid in btstrp_id_strt:btstrp_id_end
    write_dat = hcat(x_ary, LCPDF_ary_shrd[pid-btstrp_id_strt+1,:])
    write(LCPDF_file,"/$(pid)",permutedims(write_dat,(2,1)))
end

close(LCPDF_file)
