using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBSpublic
using DataFrames, GRUtils, Interpolations, RecursiveArrayTools, Dierckx, CSV
include("script_functions.jl") 

# load the simulation data 
df0 = collect_results(datadir("simulation", "SimulationIBS_B_s0var"))

# filter for certain values
rT_s = 10.0^(-0.4)
ϑ₀_s = 1
h_s = 1e-4
ϑchar_s = 6
guM_s = 0
P_s = 1
abstol_s = 1e-8
gT_range = -10:0.1:0
df = filter!([:rT, :ϑ₀, :h, :ϑchar, :guM, :P, :abstol, :gT] 
                => (rT, ϑ₀, h, ϑchar, guM, P, abstol, gT) -> 
                rT==rT_s &&
                ϑ₀==ϑ₀_s && 
                h==h_s && 
                ϑchar==ϑchar_s &&
                guM==guM_s && 
                P==P_s &&
                abstol==abstol_s && 
                gT in gT_range, 
                dropmissing(df0)
            )
# filter the IBS solutions (P=1 & gT=0)
df_ibs = filter([:guM, :gT, :s₀, :P, :h, :ϑchar, :rT] 
                => (guM, gT, s₀, P, h, ϑchar, rT) -> 
                guM==0 && 
                gT==0 && 
                s₀==0 && 
                P==1 &&
                h==h_s &&
                ϑchar==ϑchar_s &&
                rT==rT_s, 
                df0
                )

# calculate additional quantities
df = add_quant!(df)
df_ibs = add_quant!(df_ibs)

# divide the dataframe for the s0 values and sort for gT
s0var = sort(unique(df.s₀))
df_s0var = Array{DataFrame}(undef, length(s0var))
for i=1:length(s0var)
    df_s0var[i] = sort(filter([:s₀] => (s₀) -> s₀==s0var[i], df), [:gT])
end
# Plot
hold(false)
#colorscheme("solarized dark")
colorscheme("light")
plot(-df_s0var[1].gT, df_s0var[1].RS2./df_ibs.RS2)
hold(true)
for i=2:length(s0var)
    plot(-df_s0var[1].gT, df_s0var[i].RS2./df_ibs.RS2)
end
xlabel("Negative thermal gradient g_{T}")
ylabel("Resolution ratio Ξ (R_{S})")
title("Ξ (R_{S}) over g_{T} for r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑchar_s*30)°C")  
legend(string("s_{0}=",s0var[1]), string("s_{0}=",s0var[2]), string("s_{0}=",s0var[3]), string("s_{0}=",s0var[4]), string("s_{0}=",s0var[5]), string("s_{0}=",s0var[6]), string("s_{0}=",s0var[7]), location="upper right")
savefig(plotsdir("s0var","ResRatio_gT_$(ϑchar_s*30).svg"))
# export the plot data
df_export = DataFrame(gT=gT_range)
for i=1:length(s0var)
    df_export[!, "XiRS_s0=$(s0var[i])"] = df_s0var[i].RS2./df_ibs.RS2
end
CSV.write(plotsdir("s0var","ResRatio_s0var_$(ϑchar_s*30).csv"), df_export)
# export retention time difference ratio
df_export = DataFrame(gT=gT_range)
for i=1:length(s0var)
    df_export[!, "XiDthR_s0=$(s0var[i])"] = df_s0var[i].ΔþR./df_ibs.ΔþR
end
CSV.write(plotsdir("s0var","RTRatio_s0var_$(ϑchar_s*30).csv"), df_export)
# export peak width ratio
df_export = DataFrame(gT=gT_range)
for i=1:length(s0var)
    df_export[!, "XitauR_s0=$(s0var[i])"] = df_s0var[i].ðR2./df_ibs.ðR2
end
CSV.write(plotsdir("s0var","tauRatio_s0var_$(ϑchar_s*30).csv"), df_export)

# derivatives of the ratios
derivRSratio = deriv_res_ratio(df_s0var, df_ibs, s0var, gT_range)
derivRTratio = deriv_rt_ratio(df_s0var, df_ibs, s0var, gT_range)
derivtauRatio = deriv_tau_ratio(df_s0var, df_ibs, s0var, gT_range)
# export relative derivatives of the ratios 
df_rel_deriv_resRatio = DataFrame(gT=gT_range)
df_rel_deriv_RTRatio = DataFrame(gT=gT_range)
df_rel_deriv_tauRatio = DataFrame(gT=gT_range)
for i=1:length(s0var)
    df_rel_deriv_resRatio[!, "s0_$(s0var[i])"]=derivRSratio[i](gT_range)./(df_s0var[i].RS2./df_ibs.RS2)
    df_rel_deriv_RTRatio[!, "s0_$(s0var[i])"]=derivRTratio[i](gT_range)./(df_s0var[i].ΔþR./df_ibs.ΔþR)
    df_rel_deriv_tauRatio[!, "s0_$(s0var[i])"]=derivtauRatio[i](gT_range)./(df_s0var[i].ðR2./df_ibs.ðR2)
end
CSV.write(plotsdir("s0var","relDeriv_ResRatio_s0var_$(ϑchar_s*30).csv"), df_rel_deriv_resRatio)
CSV.write(plotsdir("s0var","relDeriv_RTRatio_s0var_$(ϑchar_s*30).csv"), df_rel_deriv_RTRatio)
CSV.write(plotsdir("s0var","relDeriv_tauRatio_s0var_$(ϑchar_s*30).csv"), df_rel_deriv_tauRatio)

# extimate maxima of Ξ(RS)
df_optima = optima_of_ratios(df_s0var, df_ibs, s0var)
CSV.write(plotsdir("s0var","ResRatio_s0_$(ϑchar_s*30)_optima.csv"), df_optima)
