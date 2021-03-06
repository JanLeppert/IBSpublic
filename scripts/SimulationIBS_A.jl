# Simulation of IBS for different 
# - dimensionless heating rates rT
# - dimensionless gradients gT
# - initial zone widths s₀
# - linear gradients of mobile phase velocity guM
# - gas decompression-like gradients of mobile phase velocity P
# - different solute triplets ϑchar
# Version A: pure IBS (gT=0, s₀=0, guM=0, P=1) as reference system
using DrWatson
@quickactivate "IBSpublic"
DrWatson.greet()
using IBSpublic
using DataFrames

params = Dict(
                :gT => collect(-10:0.1:0),
                :rT => 10.0^(-0.4), #.^(-3:0.1:1),
                :ϑ₀ => 1,
                :h => 1e-4,
                :ϑchar => 4,
                :char => NaN,
                :C => 0,
                :Δϑchar => 1/300,
                :s₀ => 0,
                :guM => 0,
                :P => [1, 10, 100, 1000],
                :retention_model => "ideal",        # ["ideal", "linear"]
                :retention_difference => "Δϑchar",  # ["Δϑchar", "Δlnk", "Δk"]
                :char_model => "f(ϑchar)",          # ["f(ϑchar)", "value"]
                :comparison => "ξ₀"                 # ["ξ₀", "ϑ(þ₀)", "no adjustment"]
)

Ndict = dict_list_count(params)
dicts = dict_list(params)

function makesim(d::Dict)
    @unpack gT, rT, ϑ₀, h, ϑchar, char, C, Δϑchar, s₀, guM, P, retention_model, retention_difference, char_model, comparison = d
    # initilize the parameter structure
    sys = IBSpublic.System(rT, gT, ϑ₀, guM, P)
    if char_model=="f(ϑchar)" 
        sub = IBSpublic.Substance(h, ϑchar, NaN, C, Δϑchar, s₀)
    elseif char_model=="value"
        sub = IBSpublic.Substance(h, ϑchar, char, C, Δϑchar, s₀)
    end
    options = IBSpublic.Options(retention_model, retention_difference, char_model, comparison)
    par = IBSpublic.ParTripletIBS(sys,sub,options)
    # simulation
    solþξ = Array{Any}(undef, 3)
    sols²ξ = Array{Any}(undef, 3)
    solð²ξ = Array{Any}(undef, 3)
    if comparison=="ξ₀"
        #-----------------------------------------------------------------------------------------------------------------
        # use the solution of IBS of solute "0" and adapt the temperature program of non-IBS to have the same 
        # trajectory ξ₀(þ) of the middle solute
        #----------------------------------------------------------------------------------------------------------------- 
        if gT==0 && guM==0 && P==1 && s₀==0
            for i=1:3
                solþξ[i] = IBSpublic.solving_migration(par, n=i-2, ξ₀=zero)
                sols²ξ[i] = IBSpublic.solving_bandvariance(solþξ[i], par, n=i-2, ξ₀=zero)
                solð²ξ[i] = IBSpublic.solving_peakvariance(solþξ[i], par, n=i-2, ξ₀=zero)
            end
            ξ₀func = IBSpublic.trajectory(solþξ[2], ntraj=10000)
        else
            # solve for non-IBS with adapted temperature program with
            # trajectory of solute "0" in the IBS
            # Ideal Basic Separation Solution (gT=0, s₀=0, guM=0, P=0) of the middle solute "0"
            # initilize the parameter structure for IBS
            sys_IBS = IBSpublic.System(rT, 0, ϑ₀, 0, 1)
            if char_model=="f(ϑchar)" 
                sub_IBS = IBSpublic.Substance(h, ϑchar, NaN, C, Δϑchar, 0)
            elseif char_model=="value"
                sub_IBS = IBSpublic.Substance(h, ϑchar, char, C, Δϑchar, 0)
            end
            options_IBS = IBSpublic.Options(retention_model, retention_difference, char_model, comparison)
            par_IBS = IBSpublic.ParTripletIBS(sys_IBS,sub_IBS,options_IBS)
            # simulation of IBS
            solþξ_IBS = IBSpublic.solving_migration(par_IBS, n=0, ξ₀=zero)
            # trajectory
            ξ₀func = IBSpublic.trajectory(solþξ_IBS, ntraj=10000)
            # non-IBS simulation
            for i=1:3
                solþξ[i] = IBSpublic.solving_migration(par, n=i-2, ξ₀=ξ₀func)
                sols²ξ[i] = IBSpublic.solving_bandvariance(solþξ[i], par, n=i-2, ξ₀=ξ₀func)
                solð²ξ[i] = IBSpublic.solving_peakvariance(solþξ[i], par, n=i-2, ξ₀=ξ₀func)
            end
        end
    elseif comparison=="ϑ(þ₀)"
        #------------------------------------
        # what did I meant with that?
        #------------------------------------
    elseif comparison=="no adjustment"
        #-----------------------------------------------
        # no adjustment of the temperature program
        #-----------------------------------------------
        for i=1:3
            solþξ[i] = IBSpublic.solving_migration(par, n=i-2)
            sols²ξ[i] = IBSpublic.solving_bandvariance(solþξ[i], par, n=i-2)
            solð²ξ[i] = IBSpublic.solving_peakvariance(solþξ[i], par, n=i-2)
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
    fulld[:þM] = IBSpublic.holduptime(guM, P)
    fulld[:ξ₀func] = ξ₀func
    return fulld
end

for (i, d) in enumerate(dicts)
    f = makesim(d)
    @tagsave(datadir("simulation", "SimulationIBS_A", savename(d, "bson", scientific=3, ignores=(:char, :C, :retention_model, :retention_difference, :char_model, :comparison))), f, safe=false)
end
println("script finished.")
