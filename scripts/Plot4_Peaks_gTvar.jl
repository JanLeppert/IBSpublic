using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBS
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
gT_range = -10:2:0
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
# Peaks, choosen ϑchar-values and choosen gT-values, for publication
ϑchar_s = 4
i0 = findfirst(ϑcharvar.==ϑchar_s)
gT_s = -6
j0 = findfirst(gT_range.==gT_s)
# time range 
þstart = df_ϑcharvar[i0].þRa[j0] - 6*df_ϑcharvar[i0].ðRa[j0]
þend = df_ϑcharvar[i0].þRb[j0] + 6*df_ϑcharvar[i0].ðRb[j0]
þrange = þstart:(þend-þstart)/1000:þend
# peaks
peaka = Array{Float64}(undef, length(þrange))
peakb = Array{Float64}(undef, length(þrange))
for i=1:length(þrange)
    peaka[i] = 1/(sqrt(2*π)*df_ϑcharvar[i0].ðRa[j0])*exp(-(þrange[i]-df_ϑcharvar[i0].þRa[j0])^2/(2*df_ϑcharvar[i0].ðRa[j0]^2))
    peakb[i] = 1/(sqrt(2*π)*df_ϑcharvar[i0].ðRb[j0])*exp(-(þrange[i]-df_ϑcharvar[i0].þRb[j0])^2/(2*df_ϑcharvar[i0].ðRb[j0]^2))
end
hold(false)
#colorscheme("solarized dark")
colorscheme("light")
plot(þrange, peaka)
hold(true)
plot(þrange, peakb)
xlabel("t/t_{M}")
ylabel("Abundance")
title("Peaks a and b for g_{T}=$(gT_s), r_{T}=$(round(rT_s, sigdigits=2)), T_{char,0}=$(ϑcharvar[2]*30)°C")  
legend("Peak a", "Peak b", location="lower left")
savefig(plotsdir("gTvar_rTvar","Peaks_gT=$(gT_s)_rT=$(round(rT_s, sigdigits=2))_Tchar=$(ϑcharvar[2]*30)_h=$(h_s).svg"))
# export gT and ΞRS for different gT (the data of the figure above) for one ϑchar
for i=1:length(ϑcharvar)
    df_export = DataFrame()
    for j=1:length(gT_range)
        þstart = df_ϑcharvar[i].þRa[j] - 6*df_ϑcharvar[i].ðRa[j]
        þend = df_ϑcharvar[i].þRb[j] + 6*df_ϑcharvar[i].ðRb[j]
        þrange = þstart:(þend-þstart)/1000:þend
        peaka = Array{Float64}(undef, length(þrange))
        peakb = Array{Float64}(undef, length(þrange))
        for k=1:length(þrange)
            peaka[k] = 1/(sqrt(2*π)*df_ϑcharvar[i].ðRa[j])*exp(-(þrange[k]-df_ϑcharvar[i].þRa[j])^2/(2*df_ϑcharvar[i].ðRa[j]^2))
            peakb[k] = 1/(sqrt(2*π)*df_ϑcharvar[i].ðRb[j])*exp(-(þrange[k]-df_ϑcharvar[i].þRb[j])^2/(2*df_ϑcharvar[i].ðRb[j]^2))
        end
            df_export[!, "time_gT=$(gT_range[j])"] = þrange
            df_export[!, "Peak_a_gT=$(gT_range[j])"] = peaka
            df_export[!, "Peak_b_gT=$(gT_range[j])"] = peakb
    end
    CSV.write(plotsdir("gTvar_rTvar","Peaks_Tchar=$(ϑcharvar[i]*30)_rT=$(round(rT_s, sigdigits=2))_h=$(h_s).csv"), df_export)
end
# export peaklist for different ϑchar for one gT
for i=1:length(gT_range)
    df_peaklist = peaklist(sort(filter([:gT] => (gT) -> gT==gT_range[i], df), [:ϑchar]))
    CSV.write(plotsdir("gTvar_rTvar","Peaklist_gT=$(gT_range[i])_rT=$(round(rT_s, sigdigits=2))_h=$(h_s).csv"), df_peaklist)
end