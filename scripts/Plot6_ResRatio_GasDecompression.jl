using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBS
using DataFrames, GRUtils, Interpolations, RecursiveArrayTools, Dierckx, CSV
include("script_functions.jl") 

# load the simulation data 
df0 = collect_results(datadir("simulation", "SimulationIBS_B_Pvar", "12"))

# filter for certain values
rT_s = 10.0^(-0.4)
ϑ₀_s = 1
h_s = 1e-4
ϑchar_s = 6
s₀_s = 0
guM_s = 0
#abstol_s = 1e-12
gT_range = -10:0.1:0
df = filter!([:rT, :ϑ₀, :h, :ϑchar, :s₀, :guM, :gT]#:abstol, :gT] 
                => (rT, ϑ₀, h, ϑchar, s₀, guM, gT) -> #abstol, gT) -> 
                rT==rT_s &&
                ϑ₀==ϑ₀_s && 
                h==h_s && 
                ϑchar==ϑchar_s && 
                s₀==s₀_s && 
                guM==guM_s && 
                #abstol==abstol_s && 
                gT in gT_range, 
                dropmissing(df0)
            )
# filter the IBS solutions (P=1 & gT=0)
df_ibs = filter([:guM, :gT, :s₀, :P, :ϑchar] 
                => (guM, gT, s₀, P, ϑchar) -> 
                guM==0 && 
                gT==0 && 
                s₀==0 && 
                P==1 &&
                ϑchar==ϑchar_s, 
                df0
                )

# calculate additional quantities
df = add_quant!(df)
df_ibs = add_quant!(df_ibs)

# divide the dataframe for the P values and sort for gT
Pvar = sort(unique(df.P))
df_Pvar = Array{DataFrame}(undef, length(Pvar))
for i=1:length(Pvar)
    df_Pvar[i] = sort(filter([:P] => (P) -> P==Pvar[i], df), [:gT])
end

# Plots
# Resolution, all
hold(false)
colorscheme("solarized dark")
plot(df_Pvar[1].gT, df_Pvar[1].RS2)
p_RS = gcf()
hold(true)
xlabel("Thermal Gradient g_{T}")
ylabel("Resolution R_{S}")
title("R_{S} over g_{T} for r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑchar_s*30)°C")
for i=2:length(Pvar)
    plot(df_Pvar[i].gT, df_Pvar[i].RS2)
end
draw(p_RS)
legend("1.0","1.01","1.5","2.0","2.2","3.0","5.0","100","1000","10000", location="outer upper right")
# Resolutions, choosen Pvalues
hold(false)
colorscheme("solarized dark")
plot(df_Pvar[1].gT, df_Pvar[1].RS2)
hold(true)
plot(df_Pvar[4].gT, df_Pvar[4].RS2)
plot(df_Pvar[6].gT, df_Pvar[6].RS2)
plot(df_Pvar[7].gT, df_Pvar[7].RS2)
plot(df_Pvar[8].gT, df_Pvar[8].RS2)
plot(df_Pvar[9].gT, df_Pvar[9].RS2)
xlabel("Thermal Gradient g_{T}")
ylabel("Resolution R_{S}")
title("R_{S} over g_{T} for r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑchar_s*30)°C")   
legend(string(Pvar[1]), string(Pvar[4]), string(Pvar[6]), string(Pvar[7]), string(Pvar[8]), string(Pvar[9]), location="outer upper right")
# Resolution Ratio, choosen P-values, for publication
hold(false)
#colorscheme("solarized dark")
colorscheme("light")
plot(-df_Pvar[1].gT, df_Pvar[1].RS2./df_ibs.RS2)
hold(true)
plot(-df_Pvar[4].gT, df_Pvar[4].RS2./df_ibs.RS2)
plot(-df_Pvar[6].gT, df_Pvar[6].RS2./df_ibs.RS2)
plot(-df_Pvar[7].gT, df_Pvar[7].RS2./df_ibs.RS2)
plot(-df_Pvar[8].gT, df_Pvar[8].RS2./df_ibs.RS2)
plot(-df_Pvar[9].gT, df_Pvar[9].RS2./df_ibs.RS2) 
xlabel("Negative thermal gradient -g_{T}")
ylabel("Resolution ratio Ξ (R_{S})")
title("Ξ (R_{S}) over g_{T} for r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑchar_s*30)°C")   
legend(string("P=",Pvar[1]), string("P=",Pvar[4]), string("P=",Pvar[6]), string("P=",Pvar[7]), string("P=",Pvar[8]), string("P=",Pvar[9]), location="lower left")
savefig(plotsdir("Pvar","ResRatio_GasDecompression_$(ϑchar_s*30)_full_gT.svg"))
# export gT and ΞRS for different P (the data of the figure above)
df_export = DataFrame(gT=gT_range)
for i=1:length(Pvar)
    if length(gT_range)==length(df_Pvar[i].RS2)
        df_export[!, "XiRS_P=$(Pvar[i])"] = df_Pvar[i].RS2./df_ibs.RS2
    end
end
CSV.write(plotsdir("Pvar","ResRatio_GasDecompression_$(ϑchar_s*30)_full_gT.csv"), df_export)

# plot the resolution ratio for gT=0 over Plots
df_gT0 = sort(filter([:gT] => (gT) -> gT==0, df),:P)
hold(false)
scatter(log10.(df_gT0.P), df_gT0.RS2./df_ibs.RS2)
plot(log10.(df_gT0.P), df_gT0.RS2./df_ibs.RS2)
# P=1000 ResRatio=94.28%
jG = 9/8 .* (df_gT0.P.^4 .- 1) .* (df_gT0.P.^2 .- 1) ./ (df_gT0.P.^3 .- 1).^2
scatter(jG, df_gT0.RS2./df_ibs.RS2)
scatter(log10.(df_gT0.P), jG)

# extimate maxima of Ξ(RS)
df_optima = optima_of_ratios(df_Pvar, df_ibs, Pvar)
CSV.write(plotsdir("Pvar","ResRatio_P_$(ϑchar_s*30)_optima_a.csv"), df_optima)

hold(false)
colorscheme("light")
plot(log10.(df_optima.var), df_optima.maxΞRS .- df_optima.ΞRSgT0)
plot(log10.(df_optima.var), df_optima.optgT)
plot(log10.(df_optima.var), df_optima.maxΞRS)
# export retention time difference ratio
df_export = DataFrame(gT=gT_range)
for i=1:length(Pvar)
    df_export[!, "XiDthR_P=$(Pvar[i])"] = df_Pvar[i].ΔþR./df_ibs.ΔþR
end
CSV.write(plotsdir("Pvar","RTRatio_Pvar_$(ϑchar_s*30).csv"), df_export)
# export peak width ratio
df_export = DataFrame(gT=gT_range)
for i=1:length(Pvar)
    df_export[!, "XitauR_P=$(Pvar[i])"] = df_Pvar[i].ðR2./df_ibs.ðR2
end
CSV.write(plotsdir("Pvar","tauRatio_Pvar_$(ϑchar_s*30).csv"), df_export)

# derivatives of the ratios
derivRSratio = deriv_res_ratio(df_Pvar, df_ibs, Pvar, gT_range)
derivRTratio = deriv_rt_ratio(df_Pvar, df_ibs, Pvar, gT_range)
derivtauRatio = deriv_tau_ratio(df_Pvar, df_ibs, Pvar, gT_range)
# export relative derivatives of the ratios 
df_rel_deriv_resRatio = DataFrame(gT=gT_range)
df_rel_deriv_RTRatio = DataFrame(gT=gT_range)
df_rel_deriv_tauRatio = DataFrame(gT=gT_range)
for i=1:length(Pvar)
    df_rel_deriv_resRatio[!, "P_$(Pvar[i])"]=derivRSratio[i](gT_range)./(df_Pvar[i].RS2./df_ibs.RS2)
    df_rel_deriv_RTRatio[!, "P_$(Pvar[i])"]=derivRTratio[i](gT_range)./(df_Pvar[i].ΔþR./df_ibs.ΔþR)
    df_rel_deriv_tauRatio[!, "P_$(Pvar[i])"]=derivtauRatio[i](gT_range)./(df_Pvar[i].ðR2./df_ibs.ðR2)
end
CSV.write(plotsdir("Pvar","relDeriv_ResRatio_Pvar_$(ϑchar_s*30).csv"), df_rel_deriv_resRatio)
CSV.write(plotsdir("Pvar","relDeriv_RTRatio_Pvar_$(ϑchar_s*30).csv"), df_rel_deriv_RTRatio)
CSV.write(plotsdir("Pvar","relDeriv_tauRatio_Pvar_$(ϑchar_s*30).csv"), df_rel_deriv_tauRatio)
