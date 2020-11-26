
function solving_migration(par::ParTripletIBS; n=0, ξ₀=zero, alg=OwrenZen5(), abs_tol=1e-8, rel_tol=1e-5)
    #---------------------------------------------------------------------------
    # solves the differential equations for migration þ(ξ) of a solute with the
    # parameters param of solute n and moving-gradient with trajector ξ₀(þ)
    #---------------------------------------------------------------------------
    tM = holduptime(par.sys.guM, par.sys.P)
    f_þξ(þ,p,ξ) = 1/velocity(ξ, þ, p, n=n, ξ₀=ξ₀(þ), tM=tM)
    þ₀ = 0.0
    ξspan = (0.0,1.0)
    prob_þξ = ODEProblem(f_þξ, þ₀, ξspan, par)
    solution_þξ = solve(prob_þξ, alg, abstol=abs_tol,reltol=rel_tol)
    return solution_þξ
end

function solving_bandvariance(solution_þξ, par::ParTripletIBS; n=0, ξ₀=zero, alg=OwrenZen5(), abs_tol=1e-8, rel_tol=1e-5)
    #---------------------------------------------------------------------------
    # solves the differential equations for development of the band variance
    # s²(ξ) of a solute with migration solution solution_þξ and the
    # parameters par of solute n and moving-gradient with trajector ξ₀(þ)
    # using automatic differentiation ForwardDiff
    #---------------------------------------------------------------------------
    þ(ξ) = solution_þξ(ξ)
    tM = holduptime(par.sys.guM, par.sys.P)
    υ(ξþ) = velocity(ξþ[1], ξþ[2], par, n=n, ξ₀=ξ₀(ξþ[2]), tM=tM)   # dimensionless solute velocity
    ∂υ∂ξ(ξ) = ForwardDiff.gradient(υ, [ξ,þ(ξ)])[1]      # partial deriviative for ξ of dimensionless solute velocity
    f_s²ξ(s²,p,ξ) = par.sub.h+2*s²/υ([ξ,þ(ξ)])*∂υ∂ξ(ξ)
    s²₀ = par.sub.s₀^2
    ξspan = (0.0, 1.0)
    prob_s²ξ = ODEProblem(f_s²ξ, s²₀, ξspan, par)
    solution_s²ξ = solve(prob_s²ξ, alg, abstol=abs_tol,reltol=rel_tol)
end

# add here funtion for peakvariance
function solving_peakvariance(solution_þξ, par::ParTripletIBS; n=0, ξ₀=zero, alg=OwrenZen5(), abs_tol=1e-8, rel_tol=1e-5)
    #---------------------------------------------------------------------------
    # solves the differential equations for development of the peak variance
    # ð²(ξ) of a solute with migration solution solution_þξ and the
    # parameters par of solute n and moving-gradient with trajector ξ₀(þ)
    # using automatic differentiation ForwardDiff
    #---------------------------------------------------------------------------
    þ(ξ) = solution_þξ(ξ)
    tM = holduptime(par.sys.guM, par.sys.P)
    r(ξþ) = 1/velocity(ξþ[1], ξþ[2], par, n=n, ξ₀=ξ₀(ξþ[2]), tM=tM)   # dimensionless solute residency
    ∂r∂þ(ξ) = ForwardDiff.gradient(r, [ξ,þ(ξ)])[2]      # partial deriviative for þ of dimensionless solute residency
    f_ð²ξ(ð²,p,ξ) = par.sub.h*r([ξ,þ(ξ)])^2 + 2*ð²*∂r∂þ(ξ)
    ð²₀ = par.sub.s₀^2/velocity(0, 0, par, n=n, ξ₀=ξ₀(0), tM=tM)^2
    ξspan = (0.0, 1.0)
    prob_ð²ξ = ODEProblem(f_ð²ξ, ð²₀, ξspan, par)
    solution_ð²ξ = solve(prob_ð²ξ, alg, abstol=abs_tol,reltol=rel_tol)
end

function trajectory(solution_þξ; ntraj=100, þM=1)
    # expand the solution_þξ on ntraj+1 points
    ξ = 0.0:solution_þξ.t[end]/ntraj:solution_þξ.t[end]
    þ = solution_þξ(ξ)
    # linear interpolation of the inverse solution_þξ on ntraj+1 points
    ξ₀ = LinearInterpolation((þ,), ξ, extrapolation_bc=Line())
    # if linear interpolation directly on solution_þξ, the number of points
    # can be low (<20) and the interpolation leads to differences
    return ξ₀
end