module HodgkinHuxley()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function alpha_m(V)
        return(0.1 * (V + 40.0)/(1.0 - exp(-(V + 40.0) / 10.0)))
    end

    function beta_m(V)
        return(4.0 * exp(-(V + 65.0) / 18.0))
    end

    function alpha_h(V)
        return(0.07 * exp(-(V + 65.0) / 20.0))
    end

    function beta_h(V)
        return(1.0 / (1 + exp(-(V + 35.0) / 10.0)))
    end

    function alpha_n(V)
        return(0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0)))
    end

    function beta_n(V)
        return(0.125 * exp(-(V + 65) / 80.0))
    end

    # conductances

    function gNa(g_Na, m, h)
        return(g_Na * m^(3) * h)
    end

    function gK(g_K, n)
        return(g_K*n^4)
    end

    # Channel currents

    function I_Na(V, m, h, E_Na, g_Na)
        return(gNa(g_Na, m, h)*(V - E_Na))
    end

    function I_K(V, n, E_K, g_K)
        return(gK(g_K, n) * (V - E_K))
    end

    function I_L(V, g_L, E_L)
        return(g_L * (V - E_L))
    end

    function I_stim(A_stim, I0, fStim, t)
        return(I0 .+ A_stim * sin.(2*pi*fStim*t))    
    end

    # function SingleCell!(du, u, p, t)
    #     fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
    #     du[1] = (I_stim(p[2], p[4], p[1], t) - I_Na(u[1],u[2],u[3], p[5], p[8]) - I_K(u[1],u[4], p[6], p[9]) - I_L(u[1], p[10], p[7]))/p[3]
    #     du[2] = alpha_m(u[1]) * (1.0 - u[2]) - beta_m(u[1]) * u[2]
    #     du[3] = alpha_h(u[1]) * (1.0 - u[3]) - beta_h(u[1]) * u[3]
    #     du[4] = alpha_n(u[1]) * (1.0 - u[4]) - beta_n(u[1]) * u[4]
    #     return nothing
    # end

    function SingleCell!(du, u, p, t)
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[4] .+ p[2] * sin.(2*pi*p[1]*t) - p[8]*u[2]^(3)*u[3]*(u[1] - p[5]) - p[9]*u[4]^4*(u[1] - p[6]) - p[10]*(u[1] - p[7]))/p[3]
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[3]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[3]
        du[4] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[4]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[4]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[4] .+ p[2] * sin.(2*pi*p[1]*t) - p[8]*u[2]^(3)*u[3]*(u[1] - p[5]) - p[9]*u[4]^4*(u[1] - p[6]) - p[10]*(u[1] - p[7]))/p[3]
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[3]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[3]
        du[4] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[4]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[4]
        return SA[du[1],du[2],du[3],du[4]]
    end

    function SingleCellJac!(J,u,p,t)
        J[1,1] = (-p[8] * u[3] * u[2]^3 - p[9] * u[4]^4 - p[10])/p[3]
        J[1,2] = (-3 * p[8] * u[3] * u[2]^2 * (u[1]-p[5]))/p[3]
        J[1,3] = (-p[8] * u[2]^3 * (u[1]-p[5]))/p[3]
        J[1,4] = (-4 * p[9]*u[4]^3*(u[1] - p[6]))/p[3]
        #J[1,5] = p[2] * cos.(u[5])
        J[2,1] = 0.1 * (1 - u[2])/(1 - exp(-(u[1] + 40.0) / 10.0)) + 4.0/18.0 * exp(-(u[1] + 65.0) / 18.0) * u[2] - 0.01 * exp((-40 - u[1])/10) * (1-u[2]) * (40 + u[1])/(1-exp(0.1*(-40-u[1])))^2
        J[2,2] = -4 * exp(-(u[1]+65.0)/18.0) - 0.1 * (40 + u[1])/(1 - exp(0.1 * (-40-u[1])))
        J[2,3] = 0
        J[2,4] = 0
        #J[2,5] = 0
        J[3,1] = -0.07/20.0 * exp((-65. - u[1])/20.0) * (1 - u[3]) - (0.1 * exp(0.1 * (-35 - u[1])) * u[3])/(1 + exp(0.1 * (-35.0 - u[1])))^2
        J[3,2] = 0
        J[3,3] = -0.07 * exp((-65. - u[1])/20.0) - 1/(1 + exp(0.1 * (-35 - u[1])))
        J[3,4] = 0
        #J[3,5] = 0
        J[4,1] = (0.01 * (1. - u[4]))/(1 - exp(0.1 * (-55 - u[1]))) + 0.125/80 * exp(0.0125 * (-65 - u[1])) * u[4] - (0.001 * exp(0.1 * (-55. - u[1])) * (1. - u[4]) * (55. + u[1]))/(1. - exp(0.1 * (-55. - u[1])))^2
        J[4,2] = 0
        J[4,3] = 0
        J[4,4] = -0.125 * exp(0.0125 * (-65 - u[1])) - (0.01 * (55 + u[1]))/(1 - exp(0.1 * (-55 - u[1])))
        #J[4,5] = 0     
        # J[5,1] = 0
        # J[5,2] = 0
        # J[5,3] = 0
        # J[5,4] = 0
        # J[5,5] = 0
        nothing
      end

    function Network!(du, u, p, t)
        εk, fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L, N, Vsyn, Norm, τd, τr, V0, s0, A = Tuple(p)
        for j in 1:N
            Vj, mj, hj, nj, rj = u[j:N:end]
            
            Ij = εk/Norm * (Vsyn[j] .- Vj) * sum(A[(j-1)*N+1:j*N] .* u[4*N+1:5*N])
    
            dVj = (Ij + I_stim(Astim, I0, fStim, t) - I_Na(Vj,mj,hj, E_Na, g_Na) - I_K(Vj,nj, E_K, g_K) - I_L(Vj, g_L, E_L))/Cm
            dmj = alpha_m(Vj) * (1.0 - mj) - beta_m(Vj) * mj
            dhj = alpha_h(Vj) * (1.0 - hj) - beta_h(Vj) * hj
            dnj = alpha_n(Vj) * (1.0 - nj) - beta_n(Vj) * nj
            drj = (1/τr - 1/τd)*((1 - rj)/(1+exp(-s0*(Vj - V0)))) - rj/τd
    
            du[j:N:end] = [dVj; dmj; dhj; dnj; drj]
        end    
        return nothing
    end 
end

module SynchonyMeasures

    function kuramoto_parameter_order(V, t, threshold)
        
        N = length(V[1,:])
        
        # find spiking times

        θ = zeros(length(t), N)
        raster = zeros(length(t),N)
        spks = zeros(length(t), N)
        dV = zeros(length(t)-1, N)

        for j in 1:N
            spks[:,j] = V[:,j]
            # set to zero the values bellow threshold
            for i in 1:length(t)
                if spks[i,j] > threshold
                    nothing
                else
                    spks[i,j] = 0
                end
            end

            # differential element of the potential to check for positive derivative points
            dV[:, j] = spks[2:end, j] .- spks[1:end-1, j]

            # finding the spike times as the first point where V > threshold and dV > 0 (V is increasing)
            for ti in 1:length(t)-1
                if spks[ti, j] == 0 && spks[ti+1,j] > 0
                    raster[ti,j] = 1.0
                else
                    nothing
                end
            end
            
            θj = zeros(length(t))
            indexj = findall(x -> x == 1, raster[:,j])

            for k in 1:length(indexj)-1
                for ti in indexj[k]+1:indexj[k+1]
                    θj[ti] = 2*pi*(k-1) + 2*pi*(ti - indexj[k])/(indexj[k+1]-indexj[k])
                end
            end
            θ[:,j] = θj
        end

        aux = zeros(length(θ[:,1]))
        for ti in 1:length(θ[:,1])
            aux[ti] = prod(θ[ti, :])
        end
        
        firsttime = findfirst(x -> x > 0, aux)
        lasttime = findlast(x -> x > 0, aux)

        θ = θ[firsttime:lasttime, :]
        R = zeros(length(θ[:,1]))

        for ti in 1:length(θ[:,1])
            R[ti] = abs(1/N*sum(exp.(Complex(0,1) * θ[ti,:])))
        end

        # then the kuramoto order parameter will be the mean of R:
        kuramotoOrderPar = mean(R)
        
        return(kuramotoOrderPar)
    end
end

module ISI
    using DynamicalSystems, DifferentialEquations

    function SpikeInterval(V, t, threshold)
        
        # find spiking times

        raster = zeros(length(t))
        spks = V

        # set to zero the values bellow threshold
        for i in 1:length(t)
            if spks[i] > threshold
                spks[i] = 1
            else
                spks[i] = 0
            end
        end

        # finding the spike times as the first point where V > threshold and dV > 0 (V is increasing)
        for ti in 1:length(t)-1
            if spks[ti] == 0 && spks[ti+1] > 0
                raster[ti+1] = 1.0
            else
                nothing
            end
        end

        spiketimes = t[findall(x-> x==1.0, raster)]
        isi = spiketimes[2:end] .- spiketimes[1:end-1]
        #display(plot(V[1:1000]))
        return(isi)
    end
end

module LyapunovMeasures
    import Main.HodgkinHuxley
    using DynamicalSystems, DifferentialEquations, StaticArrays, DelimitedFiles

    function rangeAfLyapunov(Amin, Amax, fMin, fMax, I0, lengthA, lengthf)
        #Constants

        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387


        # stimulation parameters
        A_stimRange = collect(Amin:round(Int64,(Amax-Amin)/lengthA):Amax)
        fStimRange = collect(fMin:(fMax-fMin)/lengthf:fMax)

        # Defining the Dynamical System

        lyaps = zeros(length(fStimRange)+1, length(A_stimRange)+1)
        
        lyaps[1, 2:end] = transpose(A_stimRange)
        lyaps[2:end, 1] = fStimRange

        # i = 1
        u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
        d0 = 1e-9
        diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-2)

        Threads.@threads for fStim in fStimRange
            println(fStim)
            Threads.@threads for Astim in A_stimRange
                A = Astim/(2 * pi * fStim * Cm)
                if A < 450
                    p0 = [fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L]
                    forcedHH = ContinuousDynamicalSystem(HodgkinHuxley.SingleCell!, u0, p0; diffeq)

                    # Simulating

                    lyap = lyapunov(forcedHH, 500; Ttr = 500, Δt = 1.0)
                    lyaps[findfirst(x -> x==fStim, fStimRange)+1,findfirst(x -> x==Astim, A_stimRange) + 1] = lyap[1]
                else
                    nothing
                end
            end
            writedlm("results/raw/lyapunov/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
        end
        return(lyaps)
    end

    function rangeAfLyapunovStat(Amin, Amax, fMin, fMax, I0, lengthA, lengthf)
         #Constants

         Cm  =   1.0 # membrane capacitance, in uF/cm^2
         g_Na = 120.0 # maximum conducances, in mS/cm^2
         g_K  =  36.0
         g_L  =   0.3
         E_Na =  50.0 # Nernst reversal potentials, in mV
         E_K  = -77.0
         E_L  = -54.387
 
 
         # stimulation parameters
         A_stimRange = collect(Amin:round(Int64,(Amax-Amin)/lengthA):Amax)
         fStimRange = collect(fMin:(fMax-fMin)/lengthf:fMax)
 
         # Defining the Dynamical System
 
         lyaps = zeros(length(fStimRange)+1, length(A_stimRange)+1)
         
         lyaps[1, 2:end] = transpose(A_stimRange)
         lyaps[2:end, 1] = fStimRange
 
         u0 = SA[-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
         d0 = 1e-9
         diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4)
 
         Threads.@threads for fStim in fStimRange
             Threads.@threads for Astim in A_stimRange
                 A = Astim/(2 * pi * fStim * Cm)
                 if A < 400
                     println("running for Amin = $(Amin), Amax = $(Amax), fStim = $(fStim), Astim = $(Astim), A = $(A)")
                     p0 = [fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L]
                     fhh = ODEFunction(HodgkinHuxley.SingleCellStat; jac = HodgkinHuxley.SingleCellJac!)
                     forcedHH = CoupledODEs(fhh, u0, p0; diffeq)
                     #forcedHH = CoupledODEs(HodgkinHuxley.SingleCellStat, u0, p0; jac=HodgkinHuxley.SingleCellJac, diffeq)
 
                     # Simulating
 
                     lyap = lyapunov(forcedHH, 1000; Ttr = 100, Δt = 0.5)
                     lyaps[findfirst(x -> x==fStim, fStimRange)+1,findfirst(x -> x==Astim, A_stimRange) + 1] = lyap
                 else
                     nothing
                 end
             end
             writedlm("results/raw/lyapunov/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
         end
         return(lyaps)
    end

    function rangeALyapunov(AMin, AMax, fStim, I0, lengthA)
        #Constants

        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387


        # stimulation parameters
        A_stimRange = collect(range(AMin,AMax, lengthA))

        # Defining the Dynamical System

        lyaps = zeros(lengthA, 2)

        lyaps[:,1] = A_stimRange

        j = 0
        u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
        diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-2)

        Threads.@threads for Astim in A_stimRange 
                j += 1
                p0 = (fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L)
                forcedHH = CoupledODEs(HodgkinHuxley.SingleCell!, u0, p0; diffeq)

                # Simulating

                lyap = lyapunov(forcedHH, 500; Ttr = 500, Δt = 1.0)
                lyaps[findfirst(x -> x==Astim, A_stimRange),2] = lyap[1]
        end
        return(lyaps)
    end
end