# functions for the notebook
θ=30
gr(html_output_format=:png);
#-------------------------------------------------------------------------------
# UI for setting of parameters
function parameters(unitopt)#, uMgrad, sinit)
    if unitopt==false #"dimensionless"
        labeltxt = "dimensionless parameters"
        Ltxt = latex("L [\\textrm{m}]")
        Lvalue = "15.0"
        Htxt = latex("H/L=1/N")
        Hvalue = "0.0001"
        tMtxt = latex("t_M [\\textrm{s}]")
        tMvalue = "40.0"
		Ptxt = latex("P")
        T0txt = latex("T_0/θ")
        T0value = "1"
        rTtxt = latex("r_T")
		rTvalue = "10.0 .^(-3:1:1)"
        gTtxt = latex("g_T")
        gTvalue = "-10:1:0"
        Tchartxt = latex("T_{char,0}/θ ")
        Tcharvalue = "[-1,0,1,2,4,6]"
        ΔCptxt = latex("ΔCp [\\textrm{J}/(\\textrm{mol K})]")
        Δtxt = latex("ΔT_{char}/θ")
        Δvalue = "1/300"
		s0txt = latex("\\sigma_{init}/L")
    elseif unitopt==true #"with units"
        labeltxt = "parameters with units"
        Ltxt = latex("L [\\textrm{m}]")
        Lvalue = "15.0"
        Htxt = latex("H [\\textrm{m}]")
        Hvalue = "0.0015"
        tMtxt = latex("t_M [\\textrm{s}]")
        tMvalue = "40.0"
		Ptxt = latex("P [-]")
        T0txt = latex("T_0 [\\degree \\textrm{C}]")
        T0value = "30"
        rTtxt = latex("R_T [\\degree \\textrm{C}/\\textrm{s}]")
		rTvalue = string("(10.0 .^(-3:1:1)).*(30.0/", tMvalue, ")")
        gTtxt = latex("G_T [\\degree \\textrm{C}/\\textrm{m}]")
        gTvalue="-20:2:0"
        Tchartxt = latex("T_{char} [\\degree \\textrm{C}]")
        Tcharvalue = "[-30,0,30,60,120,180]"
        ΔCptxt = latex("ΔCp [\\textrm{J}/(\\textrm{mol K})]")
        Δtxt = latex("ΔTchar [\\degree \\textrm{C}]")
        Δvalue = "0.1"
		s0txt = latex("\\sigma_{init} [m]")
    end
    d = OrderedDict(:label => labeltxt,
				:labelsys => "system parameters:",
                :L => textbox(hint="enter number", label=Ltxt, value=Lvalue),
                :H => textbox(hint="enter number", label=Htxt, value=Hvalue),
                :tM => textbox(hint="enter number", label=tMtxt, value=tMvalue),
				:P => textbox(hint="enter number", label=Ptxt, value="1.0"),
				:labelprog => "program parameters:",
                :T0 => textbox(hint="enter number", label=T0txt, value=T0value),
                :rT => textbox(hint="enter range", label=rTtxt, value=rTvalue),
                :gT => textbox(hint="enter range", label=gTtxt, value=gTvalue),
				:labelsub => "solute parameters:",
                :Tchar => textbox(hint="enter range", label=Tchartxt, value=Tcharvalue),
                :ΔCp => textbox(hint="enter number", label=ΔCptxt, value="0.0"),# only for test
                :Δ => textbox(hint="enter number", label=Δtxt, value=Δvalue),
				:s0 => textbox(hint="enter number", label=s0txt, value="0.0")
                )
    w = Interact.Widget{:mywidget}(d)
#if unitopt==false #"dimensionless"
#	if uMgrad==true
#		if sinit==true
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
#		elseif sinit==false
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ)))
#		end
#	elseif uMgrad==false
#		if sinit==true
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
#		elseif sinit==false
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ)))
#		end
#	end
#elseif unitopt==true #"with units"
#	@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
#	if uMgrad==true
#		if sinit==true
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
#		elseif sinit==false
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ)))
#		end
#	elseif uMgrad==false
#		if sinit==true
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
#		elseif sinit==false
#			@layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ)))
#		end
#	end
#end
    @layout! w vbox(:label, hbox(vbox(:labelsys, :L, :tM, :H, :P), vbox(:labelprog, :T0, :rT, :gT), vbox(:labelsub, :Tchar, :Δ, :s0)))
	return w
end

function uiparameter()
    uipar = @manipulate for unit in checkbox(label="quantities with SI-units")#,
							#uMgrad in checkbox(label="with gradient of mobile phase velocity"),
							#sinit in checkbox(label="with initial band width")
        w = parameters(unit)#, uMgrad, sinit)
		#@layout! w vbox(:unit, :uMgrad, :sinit, observe(_))
    end
    return uipar
end

function uiparameter_conversion(uipar)
    # convert the parameters defined in uipar=uiparameter() to the dimless or unit version
    if uipar[:unit][]==false#"dimensionless" -> convert to quantities with units
        labeltxt = "parameters with units"
        Ltxt = "L = "
        Lvalue = string(round(eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        Lunit = "\\textrm{m}"
        tMtxt = "t_M = "
        tMvalue = string(round(eval(Meta.parse(uipar[][:tM][])); sigdigits=3))
        tMunit = "\\textrm{s}"
        Htxt = "H = "
        Hvalue = string(round.(eval(Meta.parse(uipar[][:H][]))*eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        Hunit = "\\textrm{m}"
        Ptxt = "P = "
        Pvalue = string(round(eval(Meta.parse(uipar[][:P][])); sigdigits=3))
        Punit = ""
        T0txt = "T_0 = "
        T0value = string(round(eval(Meta.parse(uipar[][:T0][]))*θ; sigdigits=3))
        T0unit = "\\degree \\textrm{C}"
        rTtxt = "R_T [\\degree \\textrm{C}/\\textrm{s}] = "
        rTvalue = string(round.(eval(Meta.parse(uipar[][:rT][]))*θ/eval(Meta.parse(uipar[][:tM][])); sigdigits=3))
        gTtxt = "G_T [\\degree \\textrm{C}/\\textrm{m}] = "
        gTvalue = string(round.(eval(Meta.parse(uipar[][:gT][]))*θ/eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        Tchartxt = "T_{char,0} [\\degree \\textrm{C}] = "
        Tcharvalue = string(round.(eval(Meta.parse(uipar[][:Tchar][]))*θ; sigdigits=3))
        Δtxt = "ΔT_{char} = "
        Δvalue = string(round(eval(Meta.parse(uipar[][:Δ][]))*θ; sigdigits=3))
        Δunit = "\\degree \\textrm{C}"
        s0txt = "σ_{init} = " 
        s0value = string(round(eval(Meta.parse(uipar[][:s0][]))*eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        s0unit = "\\textrm{m}"
    elseif uipar[:unit][]==true#"with units" -> convert to dimensionless quantities
        labeltxt = "dimensionless parameters"
        Ltxt = "L = "
        Lvalue = string(round(eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        Lunit = "\\textrm{m}"
        tMtxt = "t_M = "
        tMvalue = string(round(eval(Meta.parse(uipar[][:tM][])); sigdigits=3))
        tMunit = "\\textrm{s}"
        Htxt = "H/L = 1/N = "
        Hvalue = string(round.(eval(Meta.parse(uipar[][:H][]))/eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        Hunit = ""
        Ptxt = "P = "
        Pvalue = string(round(eval(Meta.parse(uipar[][:P][])); sigdigits=3))
        Punit = ""
        T0txt = "T_0/\\theta = "
        T0value = string(round(eval(Meta.parse(uipar[][:T0][]))/θ; sigdigits=3))
        T0unit = ""
        rTtxt = "r_T = "
        rTvalue = string(round.(eval(Meta.parse(uipar[][:rT][]))*eval(Meta.parse(uipar[][:tM][])/θ); sigdigits=3))
        gTtxt = "G_T = "
        gTvalue = string(round.(eval(Meta.parse(uipar[][:gT][]))*eval(Meta.parse(uipar[][:L][])/θ); sigdigits=3))
        Tchartxt = "T_{char,0}/\\theta = "
        Tcharvalue = string(round.(eval(Meta.parse(uipar[][:Tchar][]))/θ; sigdigits=3))
        Δtxt = "ΔT_{char} = "
        Δvalue = string(round(eval(Meta.parse(uipar[][:Δ][]))/θ; sigdigits=3))
        Δunit = ""
        s0txt = "σ_{init}/L = " 
        s0value = string(round(eval(Meta.parse(uipar[][:s0][]))/eval(Meta.parse(uipar[][:L][])); sigdigits=3))
        s0unit = ""
    end
    d = OrderedDict(:label => labeltxt,
				:labelsys => "system parameters:",
                :L => latex(string(Ltxt, Lvalue, Lunit)),
                :tM => latex(string(tMtxt, tMvalue, tMunit)),
                :H => latex(string(Htxt, Hvalue, Hunit)),
				:P => latex(string(Ptxt, Pvalue, Punit)),
				:labelprog => "program parameters:",
                :T0 => latex(string(T0txt, T0value, T0unit)),
                :rT => latex(string(rTtxt, rTvalue)),
                :gT => latex(string(gTtxt, gTvalue)),
				:labelsub => "solute parameters:",
                :Tchar => latex(string(Tchartxt, Tcharvalue)),
                :Δ => latex(string(Δtxt, Δvalue, Δunit)),
				:s0 => latex(string(s0txt, s0value, s0unit))
                )
    w = Interact.Widget{:mywidget}(d)
    #@layout! w vbox(:label, :labelsys, :L, :tM, :H, :P, :labelprog, :T0, :rT, :gT, :labelsub, :Tchar, :Δ, :s0)
    @layout! w vbox(:label, vbox(:labelsys, :L, :tM, :H, :P, :labelprog, :T0, :rT, :gT, :labelsub, :Tchar, :Δ, :s0))
    return w
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# make the simulation
function makesim(d::Dict)
    @unpack gT, rT, ϑ₀, h, ϑchar, char, C, Δϑchar, s₀, guM, P, retention_model, retention_difference, char_model, comparison, abstol, reltol, L, tM = d
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
    if comparison=="ξ₀"
        #-----------------------------------------------------------------------------------------------------------------
        # use the solution of IBS of solute "0" and adapt the temperature program of non-IBS to have the same 
        # trajectory ξ₀(þ) of the middle solute
        #----------------------------------------------------------------------------------------------------------------- 
        if gT==0 && s₀==0
            for i=1:3
                solþξ[i] = IBS.solving_migration(par, n=i-2, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
                sols²ξ[i] = IBS.solving_bandvariance(solþξ[i], par, n=i-2, ξ₀=zero, abs_tol=abstol, rel_tol=reltol)
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
        end
    end
    # export the results and add them to the parameters
    fulld = copy(d)
    fulld[:solþξ_a] = DataFrame(solþξ[1])
    fulld[:solþξ_0] = DataFrame(solþξ[2])
    fulld[:solþξ_b] = DataFrame(solþξ[3])
    #fulld[:sols²ξ_a] = DataFrame(sols²ξ[1])
    #fulld[:sols²ξ_0] = DataFrame(sols²ξ[2])
    #fulld[:sols²ξ_b] = DataFrame(sols²ξ[3])
    fulld[:þRa] = solþξ[1].u[end]
    fulld[:þR0] = solþξ[2].u[end]
    fulld[:þRb] = solþξ[3].u[end]
    fulld[:sRa] = sqrt(sols²ξ[1].u[end])
    fulld[:sR0] = sqrt(sols²ξ[2].u[end])
    fulld[:sRb] = sqrt(sols²ξ[3].u[end])
    fulld[:νRa] = IBS.velocity(1, solþξ[1].u[end], par, n=-1, ξ₀=ξ₀func(solþξ[1].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:νR0] = IBS.velocity(1, solþξ[2].u[end], par, n= 0, ξ₀=ξ₀func(solþξ[2].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:νRb] = IBS.velocity(1, solþξ[3].u[end], par, n=+1, ξ₀=ξ₀func(solþξ[3].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:ðRa] = sqrt(sols²ξ[1].u[end])/IBS.velocity(1, solþξ[1].u[end], par, n=-1, ξ₀=ξ₀func(solþξ[1].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:ðR0] = sqrt(sols²ξ[2].u[end])/IBS.velocity(1, solþξ[2].u[end], par, n= 0, ξ₀=ξ₀func(solþξ[2].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:ðRb] = sqrt(sols²ξ[3].u[end])/IBS.velocity(1, solþξ[3].u[end], par, n=+1, ξ₀=ξ₀func(solþξ[3].u[end]), tM=IBS.holduptime(guM, P))
    fulld[:ΔþR] = solþξ[3].u[end]-solþξ[1].u[end]
    fulld[:ðR] = (sqrt(sols²ξ[1].u[end])/IBS.velocity(1, solþξ[1].u[end], par, n=-1, ξ₀=ξ₀func(solþξ[1].u[end]), tM=IBS.holduptime(guM, P))+sqrt(sols²ξ[3].u[end])/IBS.velocity(1, solþξ[3].u[end], par, n=+1, ξ₀=ξ₀func(solþξ[3].u[end]), tM=IBS.holduptime(guM, P)))/2
    fulld[:RS] = (solþξ[3].u[end]-solþξ[1].u[end])/(4*(sqrt(sols²ξ[1].u[end])/IBS.velocity(1, solþξ[1].u[end], par, n=-1, ξ₀=ξ₀func(solþξ[1].u[end]), tM=IBS.holduptime(guM, P))+sqrt(sols²ξ[3].u[end])/IBS.velocity(1, solþξ[3].u[end], par, n=+1, ξ₀=ξ₀func(solþξ[3].u[end]), tM=IBS.holduptime(guM, P)))/2)
    fulld[:þM] = IBS.holduptime(guM, P)
    fulld[:ξ₀func] = ξ₀func
    return fulld
end

function simulation(uipar)
    # check here if the values are in a defined range
    ΔCp = Meta.parse(uipar[][:ΔCp][])
    if uipar[:unit][]==false#"dimensionless"
		L = Meta.parse(uipar[][:L][])
		tM = Meta.parse(uipar[][:tM][])
		if Meta.parse(uipar[][:H][]) > 0
			h = Meta.parse(uipar[][:H][])
		else
			error("$(Meta.parse(uipar[][:H][])) not >0. Choose a value >0 for parameter H.")
		end
		if Meta.parse(uipar[][:T0][]) >= -8
			T0 = Meta.parse(uipar[][:T0][])
		else
			error("$(Meta.parse(uipar[][:T0][])) not >=-8. Choose a value >=-8 for parameter T0.")
		end
		if minimum(eval(Meta.parse(uipar[][:rT][]))) > 0.0001 && maximum(eval(Meta.parse(uipar[][:rT][]))) <=10
			rT = eval(Meta.parse(uipar[][:rT][]))
		else
			error("$(Meta.parse(uipar[][:rT][])) not between 0.0001 and 10. Choose a range with 0.0001<=values<=10 for parameter rT.")
		end
		if minimum(eval(Meta.parse(uipar[][:gT][]))) >= -10 && maximum(eval(Meta.parse(uipar[][:gT][]))) <=0
			gT = eval(Meta.parse(uipar[][:gT][]))
		else
			error("$(Meta.parse(uipar[][:gT][])) not between -10 and 0. Choose a range with -10<=values<=0 for parameter gT.")
		end
		# range for Tchar should depend on value of T0
		if minimum(eval(Meta.parse(uipar[][:Tchar][]))) >= -2 && maximum(eval(Meta.parse(uipar[][:Tchar][]))) <=8
			ϑchar = eval(Meta.parse(uipar[][:Tchar][]))
		else
			error("$(Meta.parse(uipar[][:Tchar][])) not between -2 and 8. Choose a range with -2<=values<=8 for parameter Tchar.")
		end
		if eval(Meta.parse(uipar[][:Δ][])) > 0
			Δ = eval(Meta.parse(uipar[][:Δ][]))
		else
			error("$(Meta.parse(uipar[][:T0][])) not >0. Choose a value >0 for parameter Δ.")
		end
        if Meta.parse(uipar[][:P][]) == 1
            P = 1
        elseif Meta.parse(uipar[][:P][]) > 1 && Meta.parse(uipar[][:P][]) <= 1000
			P = [1, Meta.parse(uipar[][:P][])]
		else
			error("$(Meta.parse(uipar[][:P][])) not between 1 and 1000. Choose a value between 1 and 1000 for parameter P. P=1 is no gradient of mobile phase velocity.")
		end
        if Meta.parse(uipar[][:s0][]) == 0
            s₀ = 0
        elseif Meta.parse(uipar[][:s0][]) > 0
			s₀ = [0, Meta.parse(uipar[][:s0][])]
		else
			error("$(Meta.parse(uipar[][:s0][])) not >0. Choose a value >0 for parameter s0.")
        end
        params = Dict(
                        :gT => collect(gT),
                        :rT => collect(rT),
                        :ϑ₀ => T0,
                        :h => h, 
                        :ϑchar => collect(ϑchar),
                        :char => NaN,
                        :C => 0,
                        :Δϑchar => Δ,
                        :s₀ => s₀,
                        :guM => 0,
                        :P => P,
                        :retention_model => "ideal",
                        :retention_difference => "Δϑchar",
                        :char_model => "f(ϑchar)",
                        :comparison => "ξ₀",
                        :abstol => 1e-8,
                        :reltol => 1e-5,
                        :L => L,
                        :tM => tM
        )
    elseif uipar[:unit][]==true#"with units"
		if Meta.parse(uipar[][:L][]) > 0
			L = Meta.parse(uipar[][:L][])
		else
			error("$(Meta.parse(uipar[][:L][])) not >0. Choose a value >0 for parameter L.")
		end
		if Meta.parse(uipar[][:tM][]) > 0
			tM = Meta.parse(uipar[][:tM][])
		else
			error("$(Meta.parse(uipar[][:tM][])) not >0. Choose a value >0 for parameter tM.")
		end
		if Meta.parse(uipar[][:H][]) > 0
			h = Meta.parse(uipar[][:H][])/L
		else
			error("$(Meta.parse(uipar[][:H][])) not >0mm. Choose a value >0mm for parameter H.")
		end
		if Meta.parse(uipar[][:T0][]) >= -240
			T0 = Meta.parse(uipar[][:T0][])/θ
		else
			error("$(Meta.parse(uipar[][:T0][])) not >=-240°C. Choose a value >=-240°C for parameter T0.")
		end
		if minimum(eval(Meta.parse(uipar[][:rT][])))*tM/θ >= 0.0001 && maximum(eval(Meta.parse(uipar[][:rT][])))*tM/θ <=10
			rT = eval(Meta.parse(uipar[][:rT][]))*tM/θ
		else
			error("$(eval(Meta.parse(uipar[][:rT][]))) not between $(0.0001*θ/tM) and $(10*θ/tM). Choose a range with $(0.0001*θ/tM)<=values<=$(10*θ/tM) for parameter rT.")
		end
		if minimum(eval(Meta.parse(uipar[][:gT][]))*L/θ) >= -10 && maximum(eval(Meta.parse(uipar[][:gT][]))*L/θ) ==0
			gT = eval(Meta.parse(uipar[][:gT][]))*L/θ
		else
			error("$(Meta.parse(uipar[][:gT][])) not between $(-10*θ/L)°C/m and 0°C/m. Choose a range with $(-10*θ/L)°C/m<=values<=0°C/m for parameter gT. The last element must be 0.")
		end
		if minimum(eval(Meta.parse(uipar[][:Tchar][]))/θ) >= -2 && maximum(eval(Meta.parse(uipar[][:Tchar][]))/θ) <=8
			ϑchar = eval(Meta.parse(uipar[][:Tchar][]))/θ
		else
			error("$(Meta.parse(uipar[][:Tchar][])) not between $(-2*θ) and $(8*θ). Choose a range with $(-2*θ)<=values<=$(8*θ) for parameter Tchar.")
		end
		if Meta.parse(uipar[][:Δ][]) > 0
			Δ = Meta.parse(uipar[][:Δ][])/θ
		else
			error("$(Meta.parse(uipar[][:T0][])) not >0. Choose a value >0 for parameter Δ.")
		end
        if Meta.parse(uipar[][:P][]) == 1
            P = 1
        elseif Meta.parse(uipar[][:P][]) > 1 && Meta.parse(uipar[][:P][]) <= 1000
			P = [1, Meta.parse(uipar[][:P][])]
		else
			error("$(Meta.parse(uipar[][:P][])) not between 1 and 1000. Choose a value between 1 and 1000 for parameter P. P=1 is no gradient of mobile phase velocity.")
		end
        if Meta.parse(uipar[][:s0][])/L==0
            s₀ = 0
        elseif Meta.parse(uipar[][:s0][])/L >= 0
            s₀ = [0, Meta.parse(uipar[][:s0][])/L]
		else
			error("$(Meta.parse(uipar[][:s0][])) not >0. Choose a value >0 for parameter s0.")
        end
        params = Dict(
                        :gT => collect(gT),#.*L/θ,
                        :rT => collect(rT),#.*tM/θ,
                        :ϑ₀ => T0,#/θ,
                        :h => h,#/L, 
                        :ϑchar => collect(ϑchar),#./θ,
                        :char => NaN,
                        :C => 0,
                        :Δϑchar => Δ,#/θ,
                        :s₀ => s₀,#/L,
                        :guM => 0,
                        :P => P,
                        :retention_model => "ideal",
                        :retention_difference => "Δϑchar",
                        :char_model => "f(ϑchar)",
                        :comparison => "ξ₀",
                        :abstol => 1e-8,
                        :reltol => 1e-5,
                        :L => L,
                        :tM => tM
        )
    end
    dicts = dict_list(params)
    Ndict = dict_list_count(params)
    simresult = Array{Any}(undef, Ndict)
    for i=1:Ndict
        simresult[i] = makesim(dicts[i])
    end
    println("calculation finished.")
    return DataFrame(simresult)
end
#-------------------------------------------------------------------------------
# extract the important parameters from simres and put them in a dict
function param_dict(simres)
    gT = sort(unique(simres.gT))
    rT = sort(unique(simres.rT))
    ϑ₀ = sort(unique(simres.ϑ₀))
    h = sort(unique(simres.h))
    ϑchar = sort(unique(simres.ϑchar))
    Δϑchar = sort(unique(simres.Δϑchar))
    s₀ = sort(unique(simres.s₀))
    P = sort(unique(simres.P))
    L = sort(unique(simres.L))
    tM = sort(unique(simres.tM))
    
    params = Dict(
                        :gT => collect(gT),
                        :rT => collect(rT),
                        :ϑ₀ => ϑ₀,
                        :h => h, 
                        :ϑchar => collect(ϑchar),
                        :Δϑchar => Δϑchar,
                        :s₀ => s₀,
                        :P => P,
                        :L => L,
                        :tM => tM
        )
    return params
end
#-------------------------------------------------------------------------------
# plot chromatogram, temperature program, performance ratios and resolution
function plot_chrom_prog(simres)
    params = param_dict(simres)
    ui0=@manipulate for unit in checkbox(label="quantities with SI-unnits")
        if unit==false#"dimensionless"
            lbl_rT = latex("r_T")
            lbl_gT = latex("-g_T")
            lbl_Tchar = latex("\\vartheta_{char}")
            tfactor = 1
            τfactor = 1
            rTfactor = 1
            gTfactor = 1
            ϑcharfactor = 1
            gT0 = 0
        elseif unit==true#"with units"
            lbl_rT = latex("R_T / (1 \\degree C/s)")
            lbl_gT = latex("-G_T \\, [\\degree C/m]")
            lbl_Tchar = latex("T_{char} [\\degree C]")
            tfactor = params[:tM][]
            τfactor = params[:tM][]
            rTfactor = θ/params[:tM][]
            gTfactor = θ/params[:L][]
            ϑcharfactor = θ
            gT0 = -0.0
        end
        rT = params[:rT].*rTfactor
        gT = params[:gT].*gTfactor
        ϑchar = params[:ϑchar].*ϑcharfactor
        ui=@manipulate for srT in slider(rT, label=lbl_rT),
                        sgT in slider(sort(-gT), label=lbl_gT, value=gT0),
                        sTchar in slider(ϑchar, label=lbl_Tchar)
            pChrom = plot_chromatogram(simres, unit, srT, -sgT, sTchar)
            pTprog = plot_program(simres, unit, srT, -sgT, sTchar)
            pRatio = plot_ratio_Δt_τ(simres, unit, srT, -sgT, sTchar)
            pRes = plot_resolution(simres, unit, srT, -sgT, sTchar)
            lay = @layout [a{0.7w} grid(3,1)] #[a{0.7w} grid(3,1)]#[b{0.5h} c]]
            pplot = plot(pChrom, pTprog, pRatio, pRes, layout=lay, size=(900,700))#plot(pChrom, pTprog, pRatio, pRes, layout=lay, size=(900,700))
            #lay = @layout [a{0.7w} grid(2,1)]
            #pplot = plot(pChrom, pTprog, pRatio, layout=lay, size=(900,700))
        end
        @layout! ui vbox(:sTchar, :srT, :sgT, observe(_), string("Figure 2. Result of the simulation of IBS."))
    end
    return ui0
end
# plot the chromatogram of solute "a" and solute "b"
function plot_chromatogram(simres, unit, srT, sgT, sTchar)
    params = param_dict(simres)
    if unit==false#"dimensionless"
		tfactor = 1
        τfactor = 1
        rTfactor = 1
        gTfactor = 1
        ϑcharfactor = 1
    elseif unit==true#"with units"
		tfactor = params[:tM][]
        τfactor = params[:tM][]
        rTfactor = θ/params[:tM][]
        gTfactor = θ/params[:L][]
        ϑcharfactor = θ
    end
    # slider values
    rT = params[:rT].*rTfactor
    gT = params[:gT].*gTfactor
    ϑchar = params[:ϑchar].*ϑcharfactor
    i = findfirst(rT.==srT[])
    j = findfirst(gT.==sgT[])
    j0 = findfirst(gT.==0)
    k = findfirst(ϑchar.==sTchar[])
    if unit==false#"dimensionless"
		titletxt = latexstring("\\textrm{Fig. 2a) Chromatogram}, r_T=$(round(rT[i]; sigdigits=2)), g_T=$(gT[j]), \\vartheta_{char}=$(ϑchar[k])")
        xlabeltxt = L"t/t_M"
        tunit = ""
        Δttxt = "\\Delta t/t_M"
        τRatxt = "\\tau_{R,a}/t_M"#"τRa/tM"
        τRbtxt = "\\tau_{R,b}/t_M"#"τRb/tM"
    elseif unit==true#"with units"
		titletxt = latexstring("\\textrm{Fig. 2a) Chromatogram}, R_T=$(round(rT[i]; sigdigits=2)) \\degree \\textrm{C}/\\textrm{s}, G_T=$(gT[j])\\degree \\textrm{C}/\\textrm{m}, T_{char}=$(ϑchar[k])\\degree \\textrm{C}")
        xlabeltxt = latexstring("t \\, [\\textrm{s}]")
        tunit = "\\textrm{s}"
        Δttxt = "\\Delta t"
        τRatxt = "\\tau_{R,a}"
        τRbtxt = "\\tau_{R,b}"
    end
    # filter simres
    df_gT0 = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                                        ->  gT==params[:gT][j0] && 
                                            rT==params[:rT][i] && 
                                            ϑchar==params[:ϑchar][k] &&
                                            P==params[:P][end] &&
                                            s₀==params[:s₀][end]
                                            , simres)  
    df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres)              
    tRa = df.þRa[]*tfactor
    tR0 = df.þR0[]*tfactor
    tRb = df.þRb[]*tfactor
    τRa = df.ðRa[]*τfactor
    τR0 = df.ðR0[]*τfactor
    τRb = df.ðRb[]*τfactor
	tstart = tRa-6*τRa
    tend = tRb+6*τRb
    tstart0 = df_gT0.þRa[]*tfactor-4*df_gT0.ðRa[]*τfactor
    tend0 = df_gT0.þRb[]*tfactor+4*df_gT0.ðRb[]*τfactor
    t = (tstart:(tend-tstart)/1000:tend)
    RS = (tRb-tRa)/(2*(τRa+τRb))
	gauss(t,tR,τR) = 1/sqrt(2*π*τR^2)*exp.(-(t.-tR).^2/(2*τR^2))
    pchrom = plot(title=titletxt, titlefontsize=14)
    xlabel!(xlabeltxt)
	xlims!(tstart0,tend0)
    plot!(pchrom, t, gauss.(t, tRa, τRa), label="a", line=2)
    plot!(pchrom, t, gauss.(t, tR0, τR0), label="0", line=2)
    plot!(pchrom, t, gauss.(t, tRb, τRb), label="b", line=2)
    # add marker lines for Δt and 2τ
    plot!(pchrom, [tRa, tRa], [0, gauss(tRa, tRa, τRa)], line=(:dot, :black), label="")
    plot!(pchrom, [tRb, tRb], [0, gauss(tRb, tRb, τRb)], line=(:dot, :black), label="")
    plot!(pchrom, [tRa, tRb], [0.9*gauss(tRb, tRb, τRb),0.9*gauss(tRb, tRb, τRb)], line=(2, :black), label="")
    plot!(pchrom, [tRa-τRa, tRa+τRa], [gauss(tRa-τRa,tRa,τRa),gauss(tRa+τRa,tRa,τRa)], line=(1, :black), label="")
    plot!(pchrom, [tRb-τRb, tRb+τRb], [gauss(tRb-τRb,tRb,τRb),gauss(tRb+τRb,tRb,τRb)], line=(1, :black), label="")
    annotate!([(tR0, 0.95*gauss(tR0, tR0, τR0), Plots.text(latexstring(Δttxt,"=$(round(tRb-tRa; sigdigits=3))",tunit), 10, :center)),
                (tR0, 0.85*gauss(tR0, tR0, τR0), Plots.text(latexstring("R_S=$(round(RS; sigdigits=3))"), 10, :center)),
                ((tRa-τRa), gauss(tRa-τRa,tRa,τRa), Plots.text(latexstring(τRatxt,"=$(round(τRa; sigdigits=3))",tunit,"\\, "), 10, :right)),
                ((tRb+τRb), gauss(tRb-τRb,tRb,τRb), Plots.text(latexstring("\\, ",τRbtxt,"=$(round(τRb; sigdigits=3))",tunit), 10, :left))
                ])

    return pchrom
end
# plot the temperature program
function plot_program(simres, unit, srT, sgT, sTchar)
    params = param_dict(simres)
    if unit==false#"dimensionless"
        ylabel!(L"T/\theta")
        xlabel!(L"t/t_M")
        Tfactor = 1
        tfactor = 1
        rTfactor = 1
        gTfactor = 1
        ϑcharfactor = 1
    elseif unit==true#"with units"
        ylabel!(latexstring("T \\, [\\degree \\textrm{C}]"))
        xlabel!(latexstring("t \\, [\\textrm{s}]"))
        Tfactor = θ
        tfactor = params[:tM][]
        rTfactor = θ/params[:tM][]
        gTfactor = θ/params[:L][]
        ϑcharfactor = θ
    end
    # slider values
    rT = params[:rT].*rTfactor
    gT = params[:gT].*gTfactor
    ϑchar = params[:ϑchar].*ϑcharfactor
    i = findfirst(rT.==srT[])
    j = findfirst(gT.==sgT[])
    k = findfirst(ϑchar.==sTchar[])
    j0 = findfirst(gT.==0)

	sys = IBS.System(params[:rT][i],params[:gT][j],params[:ϑ₀][],0,params[:P][end])
	sub = IBS.Substance(params[:h][],params[:ϑchar][k],NaN,0,params[:Δϑchar][],params[:s₀][end])
	opt = IBS.Options("ideal","Δϑchar","f(ϑchar)","ξ₀")
    p = IBS.ParTripletIBS(sys,sub,opt)
    # filter simres 
    df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres) 
    ξ₀ = df.ξ₀func[]
    þa = df.solþξ_a[].value
    ξa = df.solþξ_a[].timestamp
    þ0 = df.solþξ_0[].value
    ξ0 = df.solþξ_0[].timestamp
    þb = df.solþξ_b[].value
    ξb = df.solþξ_b[].timestamp

    Ttop = Array{Float64}(undef, length(þb))
    Tbottom = Array{Float64}(undef, length(þb))
    Ta = Array{Float64}(undef, length(þa))
	T0 = Array{Float64}(undef, length(þ0))
	Tb = Array{Float64}(undef, length(þb))
    for ii=1:length(þb)
        Ttop[ii] = IBS.temperature(0,þb[ii],p,ξ₀=ξ₀(þb[ii]))-273.15/θ
        Tbottom[ii] = IBS.temperature(1,þb[ii],p,ξ₀=ξ₀(þb[ii]))-273.15/θ
        Tb[ii] = IBS.temperature(ξb[ii],þb[ii],p,ξ₀=ξ₀(þb[ii]))-273.15/θ
    end
    for ii=1:length(þa)
		Ta[ii] = IBS.temperature(ξa[ii],þa[ii],p,ξ₀=ξ₀(þa[ii]))-273.15/θ
    end
    for ii=1:length(þ0)
        T0[ii] = IBS.temperature(ξ0[ii],þ0[ii],p,ξ₀=ξ₀(þ0[ii]))-273.15/θ
    end
    pTemp = plot(title=latexstring("\\textrm{Fig. 2b) Program}"), titlefontsize=10, legend=:outertopright)
    
	xlims!(0.0, þb[end]*tfactor)
    plot!(pTemp, þb*tfactor, Ttop*Tfactor, label=latexstring("T(0)"), line=(2,:red))
    plot!(pTemp, þb*tfactor, Tbottom*Tfactor, label=latexstring("T(L)"), line=(2,:blue))
	plot!(pTemp, þa*tfactor, Ta*Tfactor, label=latexstring("T_a"), line=(1, :solid, :black))
    plot!(pTemp, þ0*tfactor, T0*Tfactor, label=latexstring("T_0"), line=(1, :dash, :black))
	plot!(pTemp, þb*tfactor, Tb*Tfactor, label=latexstring("T_b"), line=(1, :dot, :black))

    return pTemp
end
# plot the performance ratios over gT
function plot_ratio_Δt_τ(simres, unit, srT, sgT, sTchar)
    params = param_dict(simres)
    if unit==false#"dimensionless"
        xlabeltxt = latexstring("-g_T")
        rTfactor = 1
        gTfactor = 1
        ϑcharfactor = 1
    elseif unit==true#"with units"
        xlabeltxt = latexstring("-G_T [\\degree \\textrm{C}/\\textrm{m}]")
        rTfactor = θ/params[:tM][]
        gTfactor = θ/params[:L][]
        ϑcharfactor = θ
    end
    # slider values
    rT = params[:rT].*rTfactor
    gT = params[:gT].*gTfactor
    ϑchar = params[:ϑchar].*ϑcharfactor
    i = findfirst(rT.==srT[])
    k = findfirst(ϑchar.==sTchar[])
    j0 = findfirst(gT.==0)
    j1 = findfirst(gT.==sgT[])

    ΞΔþ = Array{Float64}(undef, length(params[:gT]))
    Ξð = Array{Float64}(undef, length(params[:gT]))
    ΞRS = Array{Float64}(undef, length(params[:gT]))
    # filter simres 
    df_gT0 = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j0] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][1] &&
                    s₀==params[:s₀][1]
                    , simres) 
    for j=1:length(params[:gT])
        df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres) 
        Δþ = df.þRb[] - df.þRa[]
		Δþ_gT0 = df_gT0.þRb[] - df_gT0.þRa[]
        ΞΔþ[j] = Δþ./Δþ_gT0
		ðm = (df.ðRa[] + df.ðRb[])./2
		ðm_gT0 = (df_gT0.ðRa[] + df_gT0.ðRb[])./2
		if params[:s₀][end]==0
			Ξð[j] = ðm./ðm_gT0
		else
			Ξð[j] = (ðm./ðm_gT0)./(params[:s₀][end]*100)
		end
		ΞRS[j] = df.RS[]./df_gT0.RS[]
    end
    
    plot(legend=:outertopright, title=latexstring("\\textrm{Fig. 2c) Ratios}"), titlefontsize=16)
    plot!(-params[:gT].*gTfactor, ΞΔþ, label=latexstring("\\Xi (\\Delta t)"), line=2)
	if params[:s₀][end]==0
    	plot!(-params[:gT].*gTfactor, Ξð, label=latexstring("\\Xi (\\tau)"), line=2)
	else
		plot!(-params[:gT].*gTfactor, Ξð, label=latexstring("\\Xi (\\tau)/(100 s_0)"), line=2)
	end
	plot!(-params[:gT].*gTfactor, ΞRS, label=latexstring("\\Xi (R_S)"), line=2)
	scatter!([-params[:gT][j1]*gTfactor, -params[:gT][j1]*gTfactor, -params[:gT][j1]*gTfactor], [ΞΔþ[j1], Ξð[j1], ΞRS[j1]], label="")
	plot!([-params[:gT][j1]*gTfactor, -params[:gT][j1]*gTfactor, -params[:gT][j1]*gTfactor], sort([ΞΔþ[j1], Ξð[j1], ΞRS[j1]]), line=(1,:dot, :black), label="")
    xlims!(minimum(-params[:gT].*gTfactor), maximum(-params[:gT].*gTfactor))
    ylims!(0, maximum([ΞΔþ Ξð ΞRS]))
    xlabel!(xlabeltxt)
	ylabel!(latexstring("\\textrm{Ratio}"))
end
# plot the resolution between solute "a" and solute "b" over gT
function plot_resolution(simres, unit, srT, sgT, sTchar)
    params = param_dict(simres)
    if unit==false#"dimensionless"
        xlabeltxt = latexstring("-g_T")
        rTfactor = 1
        gTfactor = 1
        ϑcharfactor = 1
    elseif unit==true#"with units"
        xlabeltxt = latexstring("-G_T [\\degree \\textrm{C}/\\textrm{m}]")
        rTfactor = θ/params[:tM][]
        gTfactor = θ/params[:L][]
        ϑcharfactor = θ
    end
    # slider values
    rT = params[:rT].*rTfactor
    gT = params[:gT].*gTfactor
    ϑchar = params[:ϑchar].*ϑcharfactor
    i = findfirst(rT.==srT[])
    k = findfirst(ϑchar.==sTchar[])
    j1 = findfirst(gT.==sgT[])

    RS = Array{Float64}(undef, length(params[:gT]))
    for j=1:length(params[:gT])
        df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres) 
		RS[j] = df.RS[]
    end

	
	pRes = plot(title=latexstring("\\textrm{Fig. 2d) Resolution} \\, R_S"), titlefontsize=11, legend=:bottomright)
	plot!(pRes, -params[:gT].*gTfactor, RS, label="")
	scatter!(pRes, [-params[:gT][j1]*gTfactor,-params[:gT][j1]*gTfactor], [RS[j1], RS[j1]], label="")
	xlims!(minimum(-params[:gT].*gTfactor), maximum(-params[:gT].*gTfactor))
	xlabel!(xlabeltxt)
	ylabel!(latexstring("R_S"))
end
#-------------------------------------------------------------------------------
# show the peaklist
function peaklist(simres)
    params = param_dict(simres)
    ui0=@manipulate for unit in checkbox(label="quantities with SI-unnits")
        if unit==false#"dimensionless"
            lbl_rT = latex("r_T")
            lbl_gT = latex("-g_T")
            gTfactor = 1
            rTfactor = 1
            gT0 = 0
        elseif unit==true#"with units"
            lbl_rT = latex("R_T / (1 \\degree C/s)")
            lbl_gT = latex("-G_T \\, [\\degree C/m]")
            gTfactor = θ/params[:L][]
            rTfactor = θ/params[:tM][]
            gT0 = -0.0
        end
        rT = params[:rT].*rTfactor
        gT = params[:gT].*gTfactor
        ui=@manipulate for  srT in slider(rT, label=lbl_rT),
                            sgT in slider(sort(-gT), label=lbl_gT, value=gT0)
            i = findfirst(rT.==srT[])
            j = findfirst(gT.==-sgT[])
            if unit==false#"dimensionless"
                peaklist_dimless(simres, i, j)
            elseif unit==true#"with units"
                peaklist_units(simres, i, j)
            end
        end
        @layout! ui vbox(:unit, :srT, :sgT, string("Table 1. Peaklist."), observe(_))
    end
    return ui0
end
 
function peaklist_dimless(simres, i, j)
    params = param_dict(simres) 
    name = Array{Float64}(undef, length(params[:ϑchar]))
    data = Array{Float64}(undef, length(params[:ϑchar]), 9)

    for k=1:length(params[:ϑchar])
        df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres) 
        name[k] = params[:ϑchar][k]
        data[k,1] = df.þRa[]     # þR
        data[k,2] = df.þR0[]
        data[k,3] = df.þRb[]
        data[k,4] = df.ðRa[]     # \ethR
        data[k,5] = df.ðR0[]
        data[k,6] = df.ðRb[] 
        data[k,7] = df.þRb[]-df.þRa[]
        data[k,8] = (df.ðRa[]+df.ðRb[])/2
        data[k,9] = (df.þRb[]-df.þRa[])/(4*(df.ðRa[]+df.ðRb[])/2)
    end
    peaks = DataFrame(ϑ_char0 = name,
                        þ_Ra = data[:,1],
                        þ_R0 = data[:,2],
                        þ_Rb = data[:,3],
                        ð_Ra = data[:,4],
                        ð_R0 = data[:,5],
                        ð_Rb = data[:,6],
                        Δþ_R = data[:,7],
                        ð_R = data[:,8],
                        R_S = data[:,9]
                        )
    return peaks
end

function peaklist_units(simres, i, j)
    params = param_dict(simres) 
    name = Array{Float64}(undef, length(params[:ϑchar]))
    data = Array{Float64}(undef, length(params[:ϑchar]), 9)

    for k=1:length(params[:ϑchar])
        df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                ->  gT==params[:gT][j] && 
                    rT==params[:rT][i] && 
                    ϑchar==params[:ϑchar][k] &&
                    P==params[:P][end] &&
                    s₀==params[:s₀][end]
                    , simres) 
        name[k] = params[:ϑchar][k]*θ
        data[k,1] = df.þRa[]*params[:tM][]     # þR
        data[k,2] = df.þR0[]*params[:tM][] 
        data[k,3] = df.þRb[]*params[:tM][] 
        data[k,4] = df.ðRa[]*params[:tM][]      # \ethR
        data[k,5] = df.ðR0[]*params[:tM][] 
        data[k,6] = df.ðRb[]*params[:tM][]  
        data[k,7] = (df.þRb[]-df.þRa[])*params[:tM][] 
        data[k,8] = (df.ðRa[]+df.ðRb[])*params[:tM][]/2
        data[k,9] = (df.þRb[]-df.þRa[])/(4*(df.ðRa[]+df.ðRb[])/2)
    end
    peaks = DataFrame(T_char0 = name,
                        t_Ra = data[:,1],
                        t_R0 = data[:,2],
                        t_Rb = data[:,3],
                        τ_Ra = data[:,4], 
                        τ_R0 = data[:,5],
                        τ_Rb = data[:,6],
                        Δt_R = data[:,7],
                        τ_R = data[:,8],
                        R_S = data[:,9]
                        )
    return peaks
end
#-------------------------------------------------------------------------------
# plot resolution ratio
function plot_resolution_ratio(simres)
    params = param_dict(simres) 
    ui0 = @manipulate for unit in checkbox(label="quantities with SI-unnits")
        if unit==false#"dimensionless"
            lbl_rT = latex("r_T")
            lbl_gT = latex("-g_T")
            gTfactor = 1
            rTfactor = 1
        elseif unit==true#"with units"
            lbl_rT = latex("R_T / (1 \\degree C/s)")
            lbl_gT = latex("-G_T \\, [\\degree C/m]")
            gTfactor = θ/params[:L][]
            rTfactor = θ/params[:tM][]
        end
        rT = params[:rT].*rTfactor
        gT = params[:gT].*gTfactor
        ui = @manipulate for srT in slider(rT, label=lbl_rT),
                                sgT in slider(sort(-gT), label=lbl_gT)
                i1 = findfirst(rT.==srT[])
                j1 = findfirst(gT.==-sgT[])
                j0 = findfirst(gT.==0)
                if unit==false#"dimensionless"
                    titlestring = latexstring("\\textrm{Fig. 1a)}\\, r_T = ", round(params[:rT][i1]*rTfactor; sigdigits=3))
                    p1 = plot(title=titlestring,
                                xlabel=L"-g_T",
                                ylabel=L"\Xi(R_S)"
                                )
                elseif unit==true#"with units"
                    titlestring = latexstring("\\textrm{Fig. 1a)}\\, R_T = ", round(params[:rT][i1]*rTfactor; sigdigits=3), " \\degree C / s")
                    p1 = plot(title=titlestring,
                                xlabel=L"-G_T \, [\degree C/m]",
                                ylabel=L"\Xi(R_S)"
                                )
                end
                xlims!(minimum(-params[:gT].*gTfactor),maximum(-params[:gT].*gTfactor))
                ΞRS = Array{Float64}(undef, length(params[:rT]), length(params[:gT]), length(params[:ϑchar]))
                for ii=1:length(params[:rT])
                    for jj=1:length(params[:gT])
                        for kk=1:length(params[:ϑchar])
                            df_ibs = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                                            ->  gT==params[:gT][j0] && 
                                                rT==params[:rT][ii] && 
                                                ϑchar==params[:ϑchar][kk] &&
                                                P==params[:P][1] &&
                                                s₀==params[:s₀][1]
                                                , simres) 
                            df = filter([:gT, :rT, :ϑchar, :P, :s₀] => (gT, rT, ϑchar, P, s₀) 
                                            ->  gT==params[:gT][jj] && 
                                                rT==params[:rT][ii] && 
                                                ϑchar==params[:ϑchar][kk] &&
                                                P==params[:P][end] &&
                                                s₀==params[:s₀][end]
                                                , simres) 
                            ΞRS[ii,jj,kk] = df.RS[]./df_ibs.RS[]
                        end
                    end
                end
                ylims!(0.0, 1.0)#round(maximum(ΞRS), RoundUp, digits=1))
                for kk=1:length(params[:ϑchar])
                    plot!(p1, -params[:gT]*gTfactor, ΞRS[i1,:,kk],label="", line=2)
                end
                if unit==false#"dimensionless"
                    titlestring = latexstring("\\textrm{Fig. 1b)}\\, g_T = ", round(params[:gT][j1]*gTfactor; sigdigits=3))
                    p2 = plot(legend=:outertopright, title=titlestring, xlabel=L"\log(r_T)")
                    xlims!(minimum(log10.(params[:rT].*rTfactor)),maximum(log10.(params[:rT])))
                    ylims!(0.0, 1.0)#round(maximum(ΞRS), RoundUp, digits=1))
                    for kk=1:length(params[:ϑchar])
                        plot!(p2, log10.(params[:rT].*rTfactor), ΞRS[:,j1,kk],
                                line=2,
                                label=latexstring((params[:ϑchar][kk])),
                                legendtitle=latexstring("\\vartheta_{char,0}:"),
                                legendfontsize=10
                                )
                    end
                elseif unit==true#"with units"
                    titlestring = latexstring("\\textrm{Fig. 1b)}\\, G_T =", (round(params[:gT][j1]*gTfactor; sigdigits=3)), "\\degree C/m")
                    p2 = plot(legend=:outertopright, title=titlestring, xlabel=L"\log(R_T/(1\degree C/s))")
                    xlims!(minimum(log10.(params[:rT].*rTfactor)),maximum(log10.(params[:rT])))
                    ylims!(0.0, 1.0)#round(maximum(ΞRS), RoundUp, digits=1))
                    for kk=1:length(params[:ϑchar])
                        plot!(p2, log10.(params[:rT].*rTfactor), ΞRS[:,j1,kk],
                                line=2,
                                label=latexstring((params[:ϑchar][kk]*θ)),
                                legendtitle=latexstring("\\vartheta_{char,0} [\\degree C]:"),
                                legendfontsize=10
                                )
                    end
                end
                plot(p1, p2, layout=2, size = (900,500))
        end
        @layout! ui vbox(:srT, :sgT, observe(_), string("Figure 1. Resolution Ratio."))
    end
    return ui0
end
#-------------------------------------------------------------------------------