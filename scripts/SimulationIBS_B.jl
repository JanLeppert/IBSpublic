# Simulation of IBS for different 
# - dimensionless heating rates rT
# - dimensionless gradients gT
# - initial zone widths s₀
# - linear gradients of mobile phase velocity guM
# - gas decompression-like gradients of mobile phase velocity P
# - different solute triplets ϑchar
# Version B: separation without thermal gradient and ideal sample introduction (gT=0, s₀=0)
#            as reference system
using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBS
using DataFrames

params = Dict(
                :gT => collect(-4.0:0.1:0),#collect(-10:0.1:0),
                :rT => 10.0^(-0.4),#10.0.^(-3:0.1:1),
                :ϑ₀ => 1,
                :h => 1e-4,
                :ϑchar => 0,#[-1,0,1,2,4,6],#collect(-1:1:2),
                :char => NaN,
                :C => 0,
                :Δϑchar => 1/300,
                :s₀ => 0,
                :guM => 0,
                :P => 1e4,
                :retention_model => "ideal",        # ["ideal", "linear"]
                :retention_difference => "Δϑchar",  # ["Δϑchar", "Δlnk", "Δk"]
                :char_model => "f(ϑchar)",          # ["f(ϑchar)", "value"]
                :comparison => "ξ₀",                # ["ξ₀", "ϑ(þ₀)", "no adjustment"]
                :abstol => 1e-8,                    # 1e-8
                :reltol => 1e-5                     # 1e-5
)

Ndict = dict_list_count(params)
dicts = dict_list(params)

function makesim(d::Dict)
    @unpack gT, rT, ϑ₀, h, ϑchar, char, C, Δϑchar, s₀, guM, P, retention_model, retention_difference, char_model, comparison, abstol, reltol = d
    # initilize the parameter structure
    sys = IBS.System(rT, gT, ϑ₀, guM, P)
    if char_model=="f(ϑchar)" 
        sub = IBS.Substance(h, ϑchar, NaN, C, Δϑchar, s₀)
    elseif char_model=="value"
        sub = IBS.Substance(h, ϑchar, char, C, Δϑchar, s₀)
    end
    options = IBS.Options(retention_model, retention_difference, char_model, comparison)
    par = IBS.ParTripletIBS(sys,sub,options)
    # simulation
    solþξ = Array{Any}(undef, 3)
    sols²ξ = Array{Any}(undef, 3)
    solð²ξ = Array{Any}(undef, 3)
    if comparison=="ξ₀"
        #-----------------------------------------------------------------------------------------------------------------
        # use the solution of IBS of solute "0" and adapt the temperature program of non-IBS to have the same 
        # trajectory ξ₀(þ) of the middle solute
        #----------------------------------------------------------------------------------------------------------------- 
        if gT==0 && s₀==0
            for i=1:3
                solþξ[i] = IBS.solving_migration(par, n=i-2, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
                sols²ξ[i] = IBS.solving_bandvariance(solþξ[i], par, n=i-2, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
                solð²ξ[i] = IBS.solving_peakvariance(solþξ[i], par, n=i-2, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
            end
            ξ₀func = IBS.trajectory(solþξ[2], ntraj=10000)
        else
            # solve for non-IBS with adapted temperature program with
            # trajectory of solute "0" in the reference separation
            # (gT=0, s₀=0) of the middle solute "0"
            # initilize the parameter structure for reference separation
            sys_ref = IBS.System(rT, 0, ϑ₀, guM, P)
            if char_model=="f(ϑchar)" 
                sub_ref = IBS.Substance(h, ϑchar, NaN, C, Δϑchar, 0)
            elseif char_model=="value"
                sub_ref = IBS.Substance(h, ϑchar, char, C, Δϑchar, 0)
            end
            options_ref = IBS.Options(retention_model, retention_difference, char_model, comparison)
            par_ref = IBS.ParTripletIBS(sys_ref,sub_ref,options_ref)
            # 1. simulation without gradient
            solþξ_ref = IBS.solving_migration(par_ref, n=0, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
            # 2. trajectory
            ξ₀func = IBS.trajectory(solþξ_ref, ntraj=10000)
            # 3. & 4. non-IBS simulation
            for i=1:3
                solþξ[i] = IBS.solving_migration(par, n=i-2, ξ₀=ξ₀func, abs_tol=abstol, rel_tol=reltol)
                sols²ξ[i] = IBS.solving_bandvariance(solþξ[i], par, n=i-2, ξ₀=ξ₀func, abs_tol=abstol, rel_tol=reltol)
                solð²ξ[i] = IBS.solving_peakvariance(solþξ[i], par, n=i-2, ξ₀=ξ₀func, abs_tol=abstol, rel_tol=reltol)
            end
        end
    elseif comparison=="ϑ(þ₀)"
        #------------------------------------
        # placeholder
        #------------------------------------
    elseif comparison=="no adjustment"
        #-----------------------------------------------
        # no adjustment of the temperature program
        #-----------------------------------------------
        for i=1:3
            solþξ[i] = IBS.solving_migration(par, n=i-2, abs_tol=abstol, rel_tol=reltol)
            sols²ξ[i] = IBS.solving_bandvariance(solþξ[i], par, n=i-2, abs_tol=abstol, rel_tol=reltol)
            solð²ξ[i] = IBS.solving_peakvariance(solþξ[i], par, n=i-2, abs_tol=abstol, rel_tol=reltol)
        end
    end
    # export the results and add them to the parameters
    fulld = copy(d)
    fulld[:solþξ_a] = DataFrame(solþξ[1])
    fulld[:solþξ_0] = DataFrame(solþξ[2])
    fulld[:solþξ_b] = DataFrame(solþξ[3])
    fulld[:sols²ξ_a] = DataFrame(sols²ξ[1])
    fulld[:sols²ξ_0] = DataFrame(sols²ξ[2])
    fulld[:sols²ξ_b] = DataFrame(sols²ξ[3])
    fulld[:solð²ξ_a] = DataFrame(solð²ξ[1])
    fulld[:solð²ξ_0] = DataFrame(solð²ξ[2])
    fulld[:solð²ξ_b] = DataFrame(solð²ξ[3])
    fulld[:þM] = IBS.holduptime(guM, P)
    fulld[:ξ₀func] = ξ₀func
    return fulld
end

for (i, d) in enumerate(dicts)
    f = makesim(d)
    @tagsave(datadir("simulation", "SimulationIBS_B_Pvar","12", savename(d, "bson", scientific=3, ignores=(:char, :C, :retention_model, :retention_difference, :char_model, :comparison))), f, safe=false)
end
println("script finished.")
