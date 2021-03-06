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
h_s = 2e-5
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
df_ibs = filter([:guM, :gT, :s₀, :P, :h] 
                => (guM, gT, s₀, P, h) -> 
                guM==0 && 
                gT==0 && 
                s₀==0 && 
                P==1 &&
                h==h_s, 
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
# Resolution, all
hold(false)
colorscheme("solarized dark")
plot(df_ϑcharvar[1].gT, df_ϑcharvar[1].RS2)
p_RS = gcf()
hold(true)
xlabel("Thermal Gradient g_{T}")
ylabel("Resolution R_{S}")
title("R_{S} over g_{T} for r_{T}=$(round(rT_s, sigdigits=2))")
for i=2:length(ϑcharvar)
    plot(df_ϑcharvar[i].gT, df_ϑcharvar[i].RS2)
end
draw(p_RS)
legend("-30","0","30","60","120","180", location="outer upper right")
# Resolution Ratio, choosen P-values, for publication
hold(false)
#colorscheme("solarized dark")
colorscheme("light")
plot(-df_ϑcharvar[1].gT, df_ϑcharvar[1].RS2./df_ϑcharvar_ibs[1].RS2)
hold(true)
for i=2:length(ϑcharvar)
    plot(-df_ϑcharvar[i].gT, df_ϑcharvar[i].RS2./df_ϑcharvar_ibs[i].RS2)
end
xlabel("Negative thermal gradient -g_{T}")
ylabel("Resolution ratio Ξ (R_{S})")
title("Ξ (R_{S}) over g_{T} for r_{T}=$(round(rT_s, sigdigits=2))")  
legend(string("T_{char,0}=",ϑcharvar[1]*30), string("T_{char,0}=",ϑcharvar[2]*30), string("T_{char,0}=",ϑcharvar[3]*30), string("T_{char,0}=",ϑcharvar[4]*30), string("T_{char,0}=",ϑcharvar[5]*30), string("T_{char,0}=",ϑcharvar[6]*30), location="lower left")
savefig(plotsdir("gTvar_rTvar","ResRatio_gT_h=$(h_s).svg"))
# export gT and ΞRS for different P (the data of the figure above)
df_export = DataFrame(gT=gT_range)
for i=1:length(ϑcharvar)
    if length(gT_range)==length(df_ϑcharvar[i].RS2)
        df_export[!, "XiRS_Tchar=$(ϑcharvar[i]*30)"] = df_ϑcharvar[i].RS2./df_ϑcharvar_ibs[i].RS2
    end
end
CSV.write(plotsdir("gTvar_rTvar","ResRatio_gT_h=$(h_s).csv"), df_export)