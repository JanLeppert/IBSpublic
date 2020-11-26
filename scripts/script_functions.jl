# 
function add_quant!(df::DataFrame)
    df[!, :þRa] = zeros(size(df)[1])
    df[!, :þR0] = zeros(size(df)[1])
    df[!, :þRb] = zeros(size(df)[1])
    df[!, :sRa] = zeros(size(df)[1])
    df[!, :sR0] = zeros(size(df)[1])
    df[!, :sRb] = zeros(size(df)[1])
    df[!, :ðRa] = zeros(size(df)[1])
    df[!, :ðR0] = zeros(size(df)[1])
    df[!, :ðRb] = zeros(size(df)[1])
    df[!, :υRa] = zeros(size(df)[1])
    df[!, :υR0] = zeros(size(df)[1])
    df[!, :υRb] = zeros(size(df)[1])
    par = Array{Any}(undef, size(df)[1])
    for i=1:size(df)[1]
        # parameter
        sys = IBS.System(df.rT[i], df.gT[i], df.ϑ₀[i], df.guM[i], df.P[i])
        sub = IBS.Substance(df.h[i], df.ϑchar[i], df.char[i], df.C[i], df.Δϑchar[i], df.s₀[i])
        opt = IBS.Options(df.retention_model[i], df.retention_difference[i], df.char_model[i], df.comparison[i])
        par[i] = IBS.ParTripletIBS(sys,sub,opt)
        # dimensionless retention times
        df.þRa[i] = df.solþξ_a[i].value[end]
        df.þR0[i] = df.solþξ_0[i].value[end]
        df.þRb[i] = df.solþξ_b[i].value[end]
        # dimensionless band width
        df.sRa[i] = sqrt(df.sols²ξ_a[i].value[end])
        df.sR0[i] = sqrt(df.sols²ξ_0[i].value[end])
        df.sRb[i] = sqrt(df.sols²ξ_b[i].value[end])
        # dimensionless peak width
        df.ðRa[i] = sqrt(df.solð²ξ_a[i].value[end])
        df.ðR0[i] = sqrt(df.solð²ξ_0[i].value[end])
        df.ðRb[i] = sqrt(df.solð²ξ_b[i].value[end])
        # elution solute velocity
        df.υRa[i] = IBS.velocity(1, df.þRa[i], par[i], n=-1, ξ₀=df.ξ₀func[i](df.þRa[i]), tM=df.þM[i])
        df.υR0[i] = IBS.velocity(1, df.þR0[i], par[i], n= 0, ξ₀=df.ξ₀func[i](df.þR0[i]), tM=df.þM[i])
        df.υRb[i] = IBS.velocity(1, df.þRb[i], par[i], n=+1, ξ₀=df.ξ₀func[i](df.þRb[i]), tM=df.þM[i])
    end
    # retention time difference
    df[!, :ΔþR] = df.þRb .- df.þRa
    # average peak width
    df[!, :ðR1] = (df.sRa./df.υRa .+ df.sRb./df.υRb)./2
    # resolution
    df[!, :RS1] = df.ΔþR./(df.ðR1 .* 4)
    # average peak width
    df[!, :ðR2] = (df.ðRa .+ df.ðRb)./2
    # resolution
    df[!, :RS2] = df.ΔþR./(df.ðR2 .* 4)

    return df
end

# peaklist, values at z=L
function peaklist(df) 
    name = Array{Float64}(undef, size(df)[1])
    data = Array{Float64}(undef, size(df)[1], 8)
    i = 0
    for k=1:size(df)[1]
        name[k] = df.ϑchar[k]*30
        data[k,1] = df.þRa[k]     # þR
        data[k,2] = df.þR0[k]
        data[k,3] = df.þRb[k]
        data[k,4] = df.ΔþR[k]
        data[k,5] = df.ðRa[k]     # \ethR
        data[k,6] = df.ðR0[k]
        data[k,7] = df.ðRb[k] 
        data[k,8] = df.ðR1[k]   
    end
    peaks = DataFrame(name = name,
                        þRa = data[:,1],
                        þR0 = data[:,2],
                        þRb = data[:,3],
                        ΔþR = data[:,4],
                        ðRa = data[:,5],
                        ðR0 = data[:,6],
                        ðRb = data[:,7],
                        ðR = data[:,8]
                        )
    return peaks
end

# extimate maxima of Ξ(RS)
function optima_of_ratios(df, df_ibs, var)
    maxΞRS = Array{Float64}(undef, length(var))
    optgT = Array{Any}(undef, length(var))
    ΞRSgT0 = Array{Any}(undef, length(var))
    for i=1:length(var)
        #if length(gT_range)==length(df[i].RS2)
            spl = Spline1D(df[i].gT, df[i].RS2./df_ibs.RS2)
            spl_derv = Spline1D(df[i].gT, derivative(spl, df[i].gT))
            if isempty(roots(spl_derv))
                optgT[i] = 0#NaN
            else
                optgT[i] = roots(spl_derv)[end]
            end
            maxΞRS[i] = spl(optgT[i])
            ΞRSgT0[i] = spl(0)
        #else
        #    optgT[i] = NaN
        #    maxΞRS[i] = NaN
        #    ΞRSgT0[i] = NaN
        #end
    end
    df_optima = DataFrame(var=var, optgT=optgT, maxΞRS=maxΞRS, ΞRSgT0=ΞRSgT0)
    return df_optima
end

# derivatives of the resolution ratio
function deriv_res_ratio(df, df_ibs, var, gT_range)
    deriv = Array{Any}(undef, length(var))
    for i=1:length(var)
        spl = Spline1D(gT_range, df[i].RS1./df_ibs.RS1)
        d_spl(x) = derivative(spl, x)
        deriv[i] = d_spl
    end
    return deriv
end
# derivatives of the retention time difference ratio
function deriv_rt_ratio(df, df_ibs, var, gT_range)
    deriv = Array{Any}(undef, length(var))
    for i=1:length(var)
        spl = Spline1D(gT_range, df[i].ΔþR./df_ibs.ΔþR)
        d_spl(x) = derivative(spl, x)
        deriv[i] = d_spl
    end
    return deriv
end
# derivatives of the peak width ratio
function deriv_tau_ratio(df, df_ibs, var, gT_range)
    deriv = Array{Any}(undef, length(var))
    for i=1:length(var)
        spl = Spline1D(gT_range, df[i].ðR1./df_ibs.ðR1)
        d_spl(x) = derivative(spl, x)
        deriv[i] = d_spl
    end
    return deriv
end