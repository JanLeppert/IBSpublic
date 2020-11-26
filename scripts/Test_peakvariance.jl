# test of the solving function for peak variance
using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBS
using DataFrames, LambertW, Dierckx, GRUtils

# first IBS solution (gT=0, s₀=0, guM=0, P=1)
sys = IBS.System(10.0^(-0.4), 0, 1, 0, 1)
sub = IBS.Substance(1e-5, 4, NaN, 0, 1/300, 0)
options = IBS.Options("ideal", "Δϑchar", "f(ϑchar)", "ξ₀")
par_IBS = IBS.ParTripletIBS(sys,sub,options)
solþξ_IBS = Array{Any}(undef, 3)
sols²ξ_IBS = Array{Any}(undef, 3)
solð²ξ_IBS = Array{Any}(undef, 3)
for i=1:3
    solþξ_IBS[i] = IBS.solving_migration(par_IBS, n=i-2, ξ₀=zero)
    sols²ξ_IBS[i] = IBS.solving_bandvariance(solþξ_IBS[i], par_IBS, n=i-2, ξ₀=zero)
    solð²ξ_IBS[i] = IBS.solving_peakvariance(solþξ_IBS[i], par_IBS, n=i-2, ξ₀=zero)
end
ξ₀func = IBS.trajectory(solþξ_IBS[2], ntraj=10000)

þa = solþξ_IBS[1].u
ξa = solþξ_IBS[1].t
sa = sqrt.(sols²ξ_IBS[1](ξa).u)
ða = sqrt.(solð²ξ_IBS[1](ξa).u)
þ0 = solþξ_IBS[2].u
ξ0 = solþξ_IBS[2].t
s0 = sqrt.(sols²ξ_IBS[2](ξ0).u)
ð0 = sqrt.(solð²ξ_IBS[2](ξ0).u)
þb = solþξ_IBS[3].u
ξb = solþξ_IBS[3].t
sb = sqrt.(sols²ξ_IBS[3](ξb).u)
ðb = sqrt.(solð²ξ_IBS[3](ξb).u)

ua = Array{Float64}(undef, length(ξa))
for i=1:length(ξa)
    ua[i] = IBS.velocity(ξa[i], þa[i], par_IBS, n=-1)
end
u0 = Array{Float64}(undef, length(ξ0))
for i=1:length(ξ0)
    u0[i] = IBS.velocity(ξ0[i], þ0[i], par_IBS, n=0)
end
ub = Array{Float64}(undef, length(ξb))
for i=1:length(ξb)
    ub[i] = IBS.velocity(ξb[i], þb[i], par_IBS, n=1)
end

# plot ð(þ) and s(þ)/u(þ)
hold(false)
colorscheme("solarized dark")
plot(þa, ða)
hold(true)
plot(þa, sa./ua)
xlabel("þ")
ylabel("ð")
title("solute a")
legend("ð", "s/u", location="upper left")

# plot s(þ) and ð(þ)u(þ)
hold(false)
colorscheme("solarized dark")
plot(þa, sa)
hold(true)
plot(þa, ða.*ua)
xlabel("þ")
ylabel("s")
title("solute a")
legend("s", "ð*u", location="upper left")