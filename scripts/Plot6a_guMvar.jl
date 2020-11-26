using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBS
using DataFrames, GRUtils, Interpolations, RecursiveArrayTools, Dierckx, CSV
include("script_functions.jl") 

# load the simulation data 
df0 = collect_results(datadir("simulation", "SimulationIBS_B_guMvar"))

# filter for certain values
rT_s = 10.0^(-0.4)
ϑ₀_s = 1
h_s = 1e-4
ϑchar_s = 4
s₀_s = 0
P_s = 1
abstol_s = 1e-8
gT_range = -10:0.1:0
df = filter!([:rT, :ϑ₀, :h, :ϑchar, :s₀, :P, :abstol, :gT] 
                => (rT, ϑ₀, h, ϑchar, s₀, P, abstol, gT) -> 
                rT==rT_s &&
                ϑ₀==ϑ₀_s && 
                h==h_s && 
                ϑchar==ϑchar_s &&
                s₀==s₀_s && 
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

# divide the dataframe for the guM values and sort for gT
guMvar = sort(unique(df.guM))
df_guMvar = Array{DataFrame}(undef, length(guMvar))
df_guMvar_ibs = Array{DataFrame}(undef, length(guMvar))
for i=1:length(guMvar)
    df_guMvar[i] = sort(filter([:guM] => (guM) -> guM==guMvar[i], df), [:gT])
    df_guMvar_ibs[i] = sort(filter([:guM] => (guM) -> guM==guMvar[i], df_ibs), [:gT])
end

# Plot

# export the plot data

# ratios and relative derivatives of the ratios

# export ratios and relative derivatives of the ratios 