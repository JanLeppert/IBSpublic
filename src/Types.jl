# dimensionless variables
# þ = t/tM      dimensionless time
# ξ = z/L       dimensionless position
struct System
    rT::Float64         # dimless heating rate rT=RT*tM/ϴ
    gT::Float64         # dimless gradient gT=GT*L/ϴ
    ϑ₀::Float64         # dimless initial temperature ϑ₀=T₀/ϴ
    guM                 # linear mobile phase velocity gradient, dimless uM = 1 + guM * (z-1/2), if it is a tuple than (υin, υout)
    P::Float64          # gas decompression-like gradient of mobile phase velocity (use only if guM=0 and P!=0)
    #um0::Float64        # mobile phase velocity, dimless = 1
    #L::Float64          # column length, dimless = 1
end

struct Substance
    h::Float64          # dimless plate height h=H/L
    ϑchar::Float64      # dimless char. temperature ϑchar=Tchar/ϴ
    char::Float64       # dimless char. thermal const. char=ϴchar/ϴ
    C::Float64          # dimless 3rd parameter C=ΔCp/R
    Δ::Float64          # dimless retention difference [Δ=Δϑchar=ΔTchar/ϴ, Δ=Δlnk, Δ=Δk], choosen by opt.retention_difference
    s₀::Float64         # dimless starting band width s₀=σ₀/L
end

struct Options
    retention_model::String # ["ideal", "linear"]
    retention_difference::String    # ["Δϑchar", "Δlnk", "Δk"]
    char_model::String  # ["f(ϑchar)", "value"]
    comparison::String  # ["ξ₀"] #, "ϑ(þ₀)", "no adjustment"]
end

struct ParTripletIBS
    sys::System
    # sys = [γ, gT, T0, uM, L] ... system parameters
    sub::Substance
    # sub = [h, ϑchar, char, C, Δ, s₀] ... substance parameters
    opt::Options
    # opt = [retention_model, retention_difference, char_model, additional_gradient, comparison]
end

struct SolTripletIBS
    þξ::SVector{3,Any}  # solution of migration t(z)
    s²ξ::SVector{3,Any} # solution of band width development σ²(z)
    par::ParTripletIBS  # parameters
end
