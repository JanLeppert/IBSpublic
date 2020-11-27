using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBSpublic
using DataFrames, GRUtils, Interpolations, RecursiveArrayTools, Dierckx, CSV
include("script_functions.jl") 

# load the simulation data 
df0 = collect_results(datadir("simulation", "SimulationIBS_B_gTvar_rTvar"))

# filter for certain values
rT_s = 10.0^(-0.4)
ϑ₀_s = 1
h_s = 1e-4
s₀_s = 0
guM_s = 0
P_s = 1
abstol_s = 1e-8
gT_range = -10:0.1:0
df = filter!([:rT, :ϑ₀, :h, :s₀, :guM, :P, :abstol, :gT] 
                => (rT, ϑ₀, h, s₀, guM, P, abstol, gT) -> 
                rT==rT_s &&
                ϑ₀==ϑ₀_s && 
                h==h_s &&
                s₀==s₀_s && 
                guM==guM_s && 
                P==P_s &&
                abstol==abstol_s && 
                gT in gT_range, 
                dropmissing(df0)
            )
# filter the IBS solutions (P=1 & gT=0)
df_ibs = filter([:guM, :gT, :s₀, :P, :h, :rT] 
                => (guM, gT, s₀, P, h, rT) -> 
                guM==0 && 
                gT==0 && 
                s₀==0 && 
                P==1 &&
                h==h_s &&
                rT==rT_s, 
                df0
                )

# calculate additional quantities
df = add_quant!(df)
df_ibs = add_quant!(df_ibs)

# divide the dataframe for the ϑchar values and sort for gT
ϑcharvar = sort(unique(df.ϑchar))
df_ϑcharvar = Array{DataFrame}(undef, length(ϑcharvar))
df_ϑcharvar_ibs = Array{DataFrame}(undef, length(ϑcharvar))
for i=1:length(ϑcharvar)
    df_ϑcharvar[i] = sort(filter([:ϑchar] => (ϑchar) -> ϑchar==ϑcharvar[i], df), [:gT])
    df_ϑcharvar_ibs[i] = sort(filter([:ϑchar] => (ϑchar) -> ϑchar==ϑcharvar[i], df_ibs), [:gT])
end

# Plots
# Ratios, choosen ϑchar-values, for publication
ϑchar_s = 4
i0 = findfirst(ϑcharvar.==ϑchar_s)
hold(false)
#colorscheme("solarized dark")
colorscheme("light")
# resolution ratio
plot(-df_ϑcharvar[i0].gT, df_ϑcharvar[i0].RS1./df_ϑcharvar_ibs[i0].RS1)
hold(true)
# retention time difference ratio
plot(-df_ϑcharvar[i0].gT, df_ϑcharvar[i0].ΔþR./df_ϑcharvar_ibs[i0].ΔþR)
# peak width ratio
plot(-df_ϑcharvar[i0].gT, df_ϑcharvar[i0].ðR1./df_ϑcharvar_ibs[i0].ðR1)
xlabel("Negative thermal gradient -g_{T}")
ylabel("Ratio Ξ")
title("Ratios Ξ over g_{T} for r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑcharvar[2]*30)°C")  
legend("Ξ (R_{S})", "Ξ (Δþ_{R})", "Ξ (ð_{R})", location="lower left")
savefig(plotsdir("gTvar_rTvar","Ratios_gT_rT=$(round(rT_s, sigdigits=2))_Tchar=$(ϑcharvar[2]*30)_h=$(h_s).svg"))
# export gT and ΞRS for different P (the data of the figure above)
df_export = DataFrame(gT=gT_range)
for i=1:length(ϑcharvar)
    if length(gT_range)==length(df_ϑcharvar[i].RS2)
        df_export[!, "XiRS_Tchar=$(ϑcharvar[i]*30)"] = df_ϑcharvar[i].RS1./df_ϑcharvar_ibs[i].RS1
        df_export[!, "XiDeltatR_Tchar=$(ϑcharvar[i]*30)"] = df_ϑcharvar[i].ΔþR./df_ϑcharvar_ibs[i].ΔþR
        df_export[!, "XitauR_Tchar=$(ϑcharvar[i]*30)"] = df_ϑcharvar[i].ðR1./df_ϑcharvar_ibs[i].ðR1
    end
end
CSV.write(plotsdir("gTvar_rTvar","Ratios_gT_rT=$(round(rT_s, sigdigits=2))_h=$(h_s).csv"), df_export)