# Loading and evaluating thr Simulation of IBS for different 
# SimulationIBS_A.jl
# - dimensionless heating rates rT
# - dimensionless gradients gT
# - initial zone widths s₀
# - linear gradients of mobile phase velocity guM
# - gas decompression-like gradients of mobile phase velocity P
# - different solute triplets ϑchar
# Version A: pure IBS (gT=0, s₀=0, guM=0, P=0) as reference system
using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBSpublic
using DataFrames, GRUtils, Interpolations, RecursiveArrayTools
include("script_functions.jl") 

# load the simulation data 
df = collect_results!(datadir("simulation", "SimulationIBS_A"))

#------------------------------------------------------------------------------------------------------
# calculate additional quantities
df = add_quant!(df)
#---------------------------------------------------------------------------------------------------
# filter the dataframe for the guM values and sort for gT
guM_var = sort(unique(df.guM))
df_guM_var = Array{DataFrame}(undef, length(guM_var))
for i=1:length(guM_var)
    df_guM_var[i] = sort(filter([:guM, :P] => (guM, P) -> guM==guM_var[i] && P==1, df), [:gT])
end
# filter the IBS solutions
df_ibs = filter([:guM, :gT, :s₀, :P] => (guM, gT, s₀, P) -> guM==0 && gT==0 && s₀==0 && P==1, df)
#--------------------------------------------------------------------------------------------------
# Plots
# Resolution
hold(false)
colorscheme("solarized dark")
plot(df_guM_var[1].gT, df_guM_var[1].RS1)
hold(true)
p1 = gcf()
xlabel("g_T")
ylabel("R_S")
title("r_{T}=$(round(df_guM_var[1].rT[1], sigdigits=2))")
for i=2:length(guM_var)
    plot!(p1, df_guM_var[i].gT, df_guM_var[i].RS1)
end
draw(p1)
legend("0.0", "1.0", "1.6", "1.8", "1.9", location="upper left")
# Resolution Ratio
hold(false)
p2 = Figure()
plot(df_guM_var[1].gT, df_guM_var[1].RS1 ./ df_ibs.RS1)
hold(true)
xlabel("g_T")
ylabel("Xi(R_{S})")
for i=2:length(guM_var)
    plot!(p2, df_guM_var[i].gT, df_guM_var[i].RS1./df_ibs.RS1)
end
draw(p2)
legend("0.0", "1.0", "1.6", "1.8", "1.9", location="upper left")
# plot þR0 over gT
hold(false)
p3 = Figure()
plot(df_guM_var[1].gT, df_guM_var[1].þR0)
hold(true)
for i=2:length(guM_var)
    plot!(p3, df_guM_var[i].gT, df_guM_var[i].þR0)
end
draw(p3)
xlabel("g_T")
ylabel("þ_{R0}")
legend("0.0", "1.0", "1.6", "1.8", "1.9", location="lower left")

#---------------------------------------------------------------------------------------------------
# filter the dataframe for the P values and sort for gT
P_var = sort(unique(df.P))
df_P_var = Array{DataFrame}(undef, length(P_var))
for i=1:length(P_var)
    df_P_var[i] = sort(filter([:guM, :P] => (guM, P) -> guM==0 && P==P_var[i], df), [:gT])
end
# filter the IBS solutions
df_ibs = filter([:guM, :gT, :s₀, :P] => (guM, gT, s₀, P) -> guM==0 && gT==0 && s₀==0 && P==1, df)
#--------------------------------------------------------------------------------------------------
# Plots
# Resolution
hold(false)
colorscheme("solarized dark")
plot(df_P_var[1].gT, df_P_var[1].RS2)
hold(true)
p4 = gcf()
xlabel("g_T")
ylabel("R_S")
title("r_{T}=$(round(df_P_var[1].rT[1], sigdigits=2))")
for i=2:length(P_var)
    plot!(p4, df_P_var[i].gT, df_P_var[i].RS2)
end
draw(p4)
legend("1", "10", "100", "1000", location="upper left")
# Resolution Ratio
hold(false)
p5 = Figure()
plot(df_P_var[1].gT, df_P_var[1].RS1 ./ df_ibs.RS1)
hold(true)
xlabel("g_T")
ylabel("Xi(R_{S})")
for i=2:length(P_var)
    plot!(p5, df_P_var[i].gT, df_P_var[i].RS1./df_ibs.RS1)
end
draw(p5)
legend("1", "10", "100", "1000", location="upper left")
# plot þR0 over gT
hold(false)
p6 = Figure()
plot(df_P_var[1].gT, df_P_var[1].þR0)
hold(true)
for i=2:length(P_var)
    plot!(p6, df_P_var[i].gT, df_P_var[i].þR0)
end
draw(p6)
xlabel("g_T")
ylabel("þ_{R0}")
legend("1", "10", "100", "1000", location="upper left")