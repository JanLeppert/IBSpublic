# script for testing the ξ₀-function and its reference system
using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBSpublic
using DataFrames, LambertW, Dierckx, GRUtils

# first IBS solution (gT=0, s₀=0, guM=0, P=1)
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 0, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_IBS = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_IBS = Array{Any}(undef, 3)
sols²ξ_IBS = Array{Any}(undef, 3)
solð²ξ_IBS = Array{Any}(undef, 3)
for i=1:3
    solþξ_IBS[i] = IBSpublic.solving_migration(par_IBS, n=i-2, ξ₀=zero)
    sols²ξ_IBS[i] = IBSpublic.solving_bandvariance(solþξ_IBS[i], par_IBS, n=i-2, ξ₀=zero)
    solð²ξ_IBS[i] = IBSpublic.solving_peakvariance(solþξ_IBS[i], par_IBS, n=i-2, ξ₀=zero)
end
ξ₀func = IBSpublic.trajectory(solþξ_IBS[2], ntraj=10000)

# plot ξ₀func(þ₀) and ξ₀(þ₀)
þ₀ = solþξ_IBS[2].u
ξ₀ = solþξ_IBS[2].t
hold(false)
colorscheme("solarized dark")
plot(þ₀, ξ₀)
hold(true)
plot(þ₀, ξ₀func(þ₀), "--r")
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("sol", "ξ_{0}func", location="upper left")

# non-IBS solution (gT=-5), using ξ₀func
sys = IBSpublic.System(10.0^(-0.4), -5, 1, 0, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_gT = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_gT = Array{Any}(undef, 3)
sols²ξ_gT = Array{Any}(undef, 3)
solð²ξ_gT = Array{Any}(undef, 3)
for i=1:3
    solþξ_gT[i] = IBSpublic.solving_migration(par_gT, n=i-2, ξ₀=ξ₀func)
    sols²ξ_gT[i] = IBSpublic.solving_bandvariance(solþξ_IBS[i], par_gT, n=i-2, ξ₀=ξ₀func)
    solð²ξ_gT[i] = IBSpublic.solving_peakvariance(solþξ_IBS[i], par_gT, n=i-2, ξ₀=ξ₀func)
end

# plot ξ₀func(þ₀) and ξ₀(þ₀) of IBS and non-IBS (gT=-5)
þ₀_IBS = solþξ_IBS[2].u
ξ₀_IBS = solþξ_IBS[2].t
þ₀_gT = solþξ_gT[2].u
ξ₀_gT = solþξ_gT[2].t
hold(false)
colorscheme("solarized dark")
plot(þ₀_IBS, ξ₀_IBS)
hold(true)
plot(þ₀_IBS, ξ₀func(þ₀_IBS), ":*g")
plot(þ₀_gT, ξ₀_gT, ".")
plot(þ₀_gT, ξ₀func(þ₀_gT), "--r")
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("sol_IBS", "ξ_{0}func_IBS", "sol_gT", "ξ_{0}func_gT", location="upper left")

# non-IBS solution (guM=1.9), using ξ₀func
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 1.9, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_guM = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_guM = Array{Any}(undef, 3)
sols²ξ_guM = Array{Any}(undef, 3)
solð²ξ_guM = Array{Any}(undef, 3)
for i=1:3
    solþξ_guM[i] = IBSpublic.solving_migration(par_guM, n=i-2, ξ₀=ξ₀func)
    sols²ξ_guM[i] = IBSpublic.solving_bandvariance(solþξ_guM[i], par_guM, n=i-2, ξ₀=ξ₀func)
    solð²ξ_guM[i] = IBSpublic.solving_peakvariance(solþξ_guM[i], par_guM, n=i-2, ξ₀=ξ₀func)
end

# plot ξ₀func(þ₀) and ξ₀(þ₀) of IBS and non-IBS (guM=1.9)
þ₀_IBS = solþξ_IBS[2].u
ξ₀_IBS = solþξ_IBS[2].t
þ₀_guM = solþξ_guM[2].u
ξ₀_guM = solþξ_guM[2].t
hold(false)
colorscheme("solarized dark")
plot(þ₀_IBS, ξ₀_IBS)
hold(true)
plot(þ₀_IBS, ξ₀func(þ₀_IBS), ":*g")
plot(þ₀_guM, ξ₀_guM, "x")
plot(þ₀_guM, ξ₀func(þ₀_guM), "--r", linewidth=2)
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("sol_IBS", "ξ_{0}func_IBS", "sol_guM", "ξ_{0}func_guM", location="upper left")

# plot T(ξ₀(þ₀),þ₀) of IBS and non-IBS (guM=1.9)
T₀_IBS = Array{Float64}(undef, length(þ₀_IBS))
for i=1:length(þ₀_IBS)
    T₀_IBS[i] = IBSpublic.temperature(ξ₀_IBS[i], þ₀_IBS[i], par_IBS)
end
T₀_guM = Array{Float64}(undef, length(þ₀_guM))
for i=1:length(þ₀_guM)
    T₀_guM[i] = IBSpublic.temperature(ξ₀_guM[i], þ₀_guM[i], par_guM, ξ₀=ξ₀func(þ₀_guM[i]))
end
hold(false)
plot(þ₀_IBS, T₀_IBS)
hold(true)
plot(þ₀_guM, T₀_guM)
xlabel("þ_{0}")
ylabel("T_{0}")
legend("IBS", "g_{u,M}=1.9", location="upper left")

# non-IBS solution (P=100), using ξ₀func
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 0, 100)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_P = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_P = Array{Any}(undef, 3)
sols²ξ_P = Array{Any}(undef, 3)
solð²ξ_P = Array{Any}(undef, 3)
for i=1:3
    solþξ_P[i] = IBSpublic.solving_migration(par_P, n=i-2, ξ₀=ξ₀func)
    sols²ξ_P[i] = IBSpublic.solving_bandvariance(solþξ_P[i], par_P, n=i-2, ξ₀=ξ₀func)
    solð²ξ_P[i] = IBSpublic.solving_peakvariance(solþξ_P[i], par_P, n=i-2, ξ₀=ξ₀func)
end

# plot ξ₀func(þ₀) and ξ₀(þ₀) of IBS and non-IBS (P=100)
þ₀_IBS = solþξ_IBS[2].u
ξ₀_IBS = solþξ_IBS[2].t
þ₀_P = solþξ_P[2].u
ξ₀_P = solþξ_P[2].t
hold(false)
colorscheme("solarized dark")
plot(þ₀_IBS, ξ₀_IBS)
hold(true)
plot(þ₀_IBS, ξ₀func(þ₀_IBS), ":*g")
plot(þ₀_P, ξ₀_P, "x")
plot(þ₀_P, ξ₀func(þ₀_P), "--r", linewidth=2)
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("sol_IBS", "ξ_{0}func_IBS", "sol_P", "ξ_{0}func_P", location="upper left")

# -> for guM!=0 or P!=1 ξ₀(þ₀)!=ξ₀func(þ₀), which follows from changed mobile phase velocity.
# at the beginning the mobile phase velocity is lower than in the case without mobile phase
# velocity gradient, while it is higher at the end of the column.
# Also the holdup time is changed.
# Conclusion: ξ₀func estimate for gT=0, but for the value guM!=0 resp. P!=1  
#-------------------------------------------------------------------------------------------------------
# non-IBS solution (guM=1.9), using ξ₀func at gT=0 and guM=1.9
# first IBS solution (gT=0, s₀=0)
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 1.9, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_gT0 = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_gT0 = Array{Any}(undef, 3)
sols²ξ_gT0 = Array{Any}(undef, 3)
solð²ξ_gT0 = Array{Any}(undef, 3)
for i=1:3
    solþξ_gT0[i] = IBSpublic.solving_migration(par_gT0, n=i-2, ξ₀=zero)
    sols²ξ_gT0[i] = IBSpublic.solving_bandvariance(solþξ_gT0[i], par_gT0, n=i-2, ξ₀=zero)
    solð²ξ_gT0[i] = IBSpublic.solving_peakvariance(solþξ_gT0[i], par_gT0, n=i-2, ξ₀=zero)
end
ξ₀func_gT0 = IBSpublic.trajectory(solþξ_gT0[2], ntraj=10000)
# solution for gT!=0
sys = IBSpublic.System(10.0^(-0.4), -5, 1, 1.9, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_guMa = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_guMa = Array{Any}(undef, 3)
sols²ξ_guMa = Array{Any}(undef, 3)
solð²ξ_guMa = Array{Any}(undef, 3)
for i=1:3
    solþξ_guMa[i] = IBSpublic.solving_migration(par_guMa, n=i-2, ξ₀=ξ₀func_gT0)
    sols²ξ_guMa[i] = IBSpublic.solving_bandvariance(solþξ_guMa[i], par_guMa, n=i-2, ξ₀=ξ₀func_gT0)
    solð²ξ_guMa[i] = IBSpublic.solving_peakvariance(solþξ_guMa[i], par_guMa, n=i-2, ξ₀=ξ₀func_gT0)
end
# plot ξ₀func(þ₀) and ξ₀(þ₀) of IBS and non-IBS (guM=1.9)
þ₀_gT0 = solþξ_gT0[2].u
ξ₀_gT0 = solþξ_gT0[2].t
þ₀_guM = solþξ_guMa[2].u
ξ₀_guM = solþξ_guMa[2].t
hold(false)
colorscheme("solarized dark")
plot(þ₀_gT0, ξ₀_gT0)
hold(true)
plot(þ₀_gT0, ξ₀func_gT0(þ₀_gT0), ":*g")
plot(þ₀_guM, ξ₀_guM, "x")
plot(þ₀_guM, ξ₀func_gT0(þ₀_guM), "--r", linewidth=2)
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("sol_IBS", "ξ_{0}func_IBS", "sol_guM", "ξ_{0}func_guM", location="upper left")

#-------------------------------------------------------------------------------------------------
# it follows from this approach, that ξ₀func for different values of guM resp. P are different:
# guM = 1.9
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 1.9, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_gT0_19 = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_gT0_19 = Array{Any}(undef, 3)
sols²ξ_gT0_19 = Array{Any}(undef, 3)
solð²ξ_gT0_19 = Array{Any}(undef, 3)
for i=1:3
    solþξ_gT0_19[i] = IBSpublic.solving_migration(par_gT0_19, n=i-2, ξ₀=zero)
    sols²ξ_gT0_19[i] = IBSpublic.solving_bandvariance(solþξ_gT0_19[i], par_gT0_19, n=i-2, ξ₀=zero)
    solð²ξ_gT0_19[i] = IBSpublic.solving_peakvariance(solþξ_gT0_19[i], par_gT0_19, n=i-2, ξ₀=zero)
end
ξ₀func_gT0_19 = IBSpublic.trajectory(solþξ_gT0_19[2], ntraj=10000)
# guM = 1.0
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 1.0, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_gT0_10 = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_gT0_10 = Array{Any}(undef, 3)
sols²ξ_gT0_10 = Array{Any}(undef, 3)
solð²ξ_gT0_10 = Array{Any}(undef, 3)
for i=1:3
    solþξ_gT0_10[i] = IBSpublic.solving_migration(par_gT0_10, n=i-2, ξ₀=zero)
    sols²ξ_gT0_10[i] = IBSpublic.solving_bandvariance(solþξ_gT0_10[i], par_gT0_10, n=i-2, ξ₀=zero)
    solð²ξ_gT0_10[i] = IBSpublic.solving_peakvariance(solþξ_gT0_10[i], par_gT0_10, n=i-2, ξ₀=zero)
end
ξ₀func_gT0_10 = IBSpublic.trajectory(solþξ_gT0_10[2], ntraj=10000)
# guM = 0
sys = IBSpublic.System(10.0^(-0.4), 0, 1, 0, 1)
sub = IBSpublic.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBSpublic.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_gT0_0 = IBSpublic.ParTripletIBS(sys,sub,options)
solþξ_gT0_0 = Array{Any}(undef, 3)
sols²ξ_gT0_0 = Array{Any}(undef, 3)
solð²ξ_gT0_0 = Array{Any}(undef, 3)
for i=1:3
    solþξ_gT0_0[i] = IBSpublic.solving_migration(par_gT0_0, n=i-2, ξ₀=zero)
    sols²ξ_gT0_0[i] = IBSpublic.solving_bandvariance(solþξ_gT0_0[i], par_gT0_0, n=i-2, ξ₀=zero)
    solð²ξ_gT0_0[i] = IBSpublic.solving_peakvariance(solþξ_gT0_0[i], par_gT0_0, n=i-2, ξ₀=zero)
end
ξ₀func_gT0_0 = IBSpublic.trajectory(solþξ_gT0_0[2], ntraj=10000)
# plot the ξ₀func
þ₀_gT0_0 = solþξ_gT0_0[2].u
þ₀_gT0_10 = solþξ_gT0_10[2].u
þ₀_gT0_19 = solþξ_gT0_19[2].u
hold(false)
plot(þ₀_gT0_0, ξ₀func_gT0_0(þ₀_gT0_0))
hold(true)
plot(þ₀_gT0_10, ξ₀func_gT0_10(þ₀_gT0_10))
plot(þ₀_gT0_19, ξ₀func_gT0_19(þ₀_gT0_19))
xlabel("þ_{0}")
ylabel("ξ_{0}")
legend("ξ_{0}func g_{u,M}=0", "ξ_{0}func g_{u,M}=1.0", "ξ_{0}func g_{u,M}=1.9", location="upper left")
# plot T₀
þ₀_gT0_0 = solþξ_gT0_0[2].u
þ₀_gT0_10 = solþξ_gT0_10[2].u
þ₀_gT0_19 = solþξ_gT0_19[2].u
ξ₀_gT0_0 = solþξ_gT0_0[2].t
ξ₀_gT0_10 = solþξ_gT0_10[2].t
ξ₀_gT0_19 = solþξ_gT0_19[2].t
T₀_gT0_0 = Array{Float64}(undef, length(þ₀_gT0_0))
for i=1:length(þ₀_gT0_0)
    T₀_gT0_0[i] = IBSpublic.temperature(ξ₀_gT0_0[i], þ₀_gT0_0[i], par_gT0_0, ξ₀=ξ₀func_gT0_0(þ₀_gT0_0[i]))
end
T₀_gT0_10 = Array{Float64}(undef, length(þ₀_gT0_10))
for i=1:length(þ₀_gT0_10)
    T₀_gT0_10[i] = IBSpublic.temperature(ξ₀_gT0_10[i], þ₀_gT0_10[i], par_gT0_10, ξ₀=ξ₀func_gT0_10(þ₀_gT0_10[i]))
end
T₀_gT0_19 = Array{Float64}(undef, length(þ₀_gT0_19))
for i=1:length(þ₀_gT0_19)
    T₀_gT0_19[i] = IBSpublic.temperature(ξ₀_gT0_19[i], þ₀_gT0_19[i], par_gT0_19, ξ₀=ξ₀func_gT0_19(þ₀_gT0_19[i]))
end
hold(false)
plot(þ₀_gT0_0, T₀_gT0_0)
hold(true)
plot(þ₀_gT0_10, T₀_gT0_10)
plot(þ₀_gT0_19, T₀_gT0_19)
xlabel("þ_{0}")
ylabel("T_{0}")
legend("g_{u,M}=0", "g_{u,M}=1.0", "g_{u,M}=1.9", location="upper left")