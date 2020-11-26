const θb = 30.0
const ϑst = 273.15/θb

#-------------------------------------------------------------------------------
function temperature(ξ, þ, p::ParTripletIBS; ξ₀=0)
        # for moving gradient input for z₀ following z₀(t) the trajector of the
        # moving gradient
        ϑ = p.sys.ϑ₀ + p.sys.gT*(ξ-ξ₀) + p.sys.rT*þ + ϑst
        return ϑ
end

function retention_factor(ξ, þ, p::ParTripletIBS; n=0, ξ₀=0)
    ϑ = temperature(ξ, þ, p, ξ₀=ξ₀)
    ϕ = 0.001 # dimensionless film thickness df/d, add ϕ to the system parameters?
    θst = 22*(1000*ϕ)^0.09
    msg_retention_model = "No such retention_model option. Choose 'ideal' or 'linear'."
    msg_retention_difference = "No such retention_difference option. Choose 'Δϑchar', 'Δlnk' or 'Δk'."
    msg_char_model = "No such char_model option. Choose 'f(ϑchar)' or 'value'."
    if p.opt.retention_difference=="Δϑchar" # Δ=ΔTchar, 3/2-param-model
        ϑchar = p.sub.ϑchar + n*p.sub.Δ + ϑst
        if p.opt.char_model=="f(ϑchar)"
            char = θst/θb*(ϑchar/ϑst)^0.7
        elseif p.opt.char_model=="value"
            char = p.sub.char
        else
            println(msg_char_model)
            return
        end
        if p.opt.retention_model=="ideal"
            lnk = (p.sub.C + ϑchar/char) * (ϑchar/ϑ - 1) + p.sub.C*log(ϑ/ϑchar)
        elseif p.opt.retention_model=="linear"
            lnk = -(ϑ - ϑchar)/char
        else
            println(msg_retention_model)
            return
        end
        k = exp(lnk)
    elseif p.opt.retention_difference=="Δlnk" # Δ=Δlnk
        ϑchar = p.sub.ϑchar + ϑst
        if p.opt.char_model=="f(ϑchar)"
            char = θst/θb*(ϑchar/ϑst)^0.7
        elseif p.opt.char_model=="value"
            char = p.sub.char
        else
            println(msg_char_model)
            return
        end
        if p.opt.retention_model=="ideal"
            lnk = (p.sub.C + ϑchar/char) * (ϑchar/ϑ - 1) + p.sub.C*log(ϑ/ϑchar) + n*p.sub.Δ
        elseif p.opt.retention_model=="linear"
            lnk = -(ϑ-ϑchar)/char + n*p.sub.Δ
        else
            println(msg_retention_model)
            return
        end
        k = exp(lnk)
    elseif p.opt.retention_difference=="Δk"
        ϑchar = p.sub.ϑchar + ϑst
        if p.opt.char_model=="f(ϑchar)"
            char = θst/θb*(ϑchar/ϑst)^0.7
        elseif p.opt.char_model=="value"
            char = p.sub.char
        else
            println(msg_char_model)
            return
        end
        if p.opt.retention_model=="ideal"
            lnk = (p.sub.C + ϑchar/char) * (ϑchar/ϑ - 1) + p.sub.C*log(ϑ/ϑchar)
        elseif p.opt.retention_model=="linear"
            lnk = -(ϑ - ϑchar)/char
        else
            println(msg_retention_model)
            return
        end
        k = exp(lnk) + n*p.sub.Δ
    else
        println(msg_retention_difference)
        return
    end
    return k
end

function mobility(ξ, þ, p::ParTripletIBS; n=0, ξ₀=0)
    # μ = u/uM, mobility is identical to the dimensionless solute velocity
    # if no gradient of the mobile phase velocity exists
    μ = 1/(1 + retention_factor(ξ, þ, p, n=n, ξ₀=ξ₀))
    return μ
end

function holduptime(guM, P)
    # dimensionless hold-up-time þM 
    # þM = tM/(L/uM)
    if P==1
        # linear mobile phase velocity gradient
        if isa(guM, Tuple)==true
            þM = 1
        else
            if guM!=0
                þM = (log(1+guM/2)-log(1-guM/2))/guM
            elseif guM==0
                þM = 1
            end
        end
    elseif P!=1
        # gas decompression-like mobile phase velocity gradient 
        # James-Martin-Compressibility
        jM = 3/2*(P^2-1)/(P^3-1)
        þM = 1/jM
    end
    return þM
end

function mobile_phase_velocity(ξ, p::ParTripletIBS; tM=1)
    if p.sys.P==1
        # dimless mobile phase velocity with a linear gradient
        # inverse of tm is the average (temporal) mobile phase velocity
        if isa(p.sys.guM, Tuple)==true
            υin = p.sys.guM[1]
            υout = p.sys.guM[2]
            υM = υin + (υout-υin)*ξ
        else
            if p.sys.guM==0
                υM = 1
            else
                υM = (1 + p.sys.guM*(ξ-1/2))*tM
            end
        end
    else
        # gas decompression-like mobile phase velocity gradient
        # James-Martin-Compressibility
        jM = 3/2*(p.sys.P^2-1)/(p.sys.P^3-1)
        υM = 1/(jM*sqrt(p.sys.P^2-ξ*(p.sys.P^2-1)))
    end
    return υM #\upsilon m
end

function velocity(ξ, þ, p::ParTripletIBS; n=0, ξ₀=0, tM=1)
    υ = mobility(ξ, þ, p, n=n, ξ₀=ξ₀)*mobile_phase_velocity(ξ, p, tM=tM)
    return υ #\upsilon
end
