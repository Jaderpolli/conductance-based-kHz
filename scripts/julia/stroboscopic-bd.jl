include("../../fhh/julia/forcedHH.jl")

import Main.HodgkinHuxley
import Main.ISI

using DelimitedFiles, DifferentialEquations, DynamicalSystems, BenchmarkTools

function bifurcationdiag()
    #Constants

    Cm  =   1.0 # membrane capacitance, in uF/cm^2
    g_Na = 120.0 # maximum conducances, in mS/cm^2
    g_K  =  36.0
    g_L  =   0.3
    E_Na =  50.0 # Nernst reversal potentials, in mV
    E_K  = -77.0
    E_L  = -54.387


    # stimulation parameters
    lengthA = 2000
    #fStim = 1.65
    I0Values = [20]
    fStimValues = [5.0]
    AminValues = [10]
    AmaxValues = [1000]

    u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
    diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-2)
    T = 1000
    thresholdm3h = 0.12
    mkpath("results/raw/bd")
    mkpath("results/raw/isi")

    for I0 in I0Values
        j = 0
        for fStim in fStimValues
            j += 1
            Amin = AminValues[j]
            Amax = AmaxValues[j]
            A_stimRange = collect(Amin:(Amax-Amin)/lengthA:Amax)
            size = round(Int64, T*fStim)
            bd = hcat(A_stimRange, zeros(length(A_stimRange), size+1))
            bdm3h = Array{Union{Nothing,Matrix{Float64}}}(nothing,length(A_stimRange))
            Threads.@threads for k in eachindex(A_stimRange)
                Astim = A_stimRange[k]
                p0 = [fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L]
                forcedHH = ContinuousDynamicalSystem(HodgkinHuxley.SingleCell!, u0, p0; diffeq)
                traj, t = trajectory(forcedHH, T; Δt = 1/fStim, Ttr = 2000)
                bd[k, 2:end] = transpose(Matrix(traj)[:,1])
        
                # Y, t = trajectory(forcedHH, T; Ttr = 1000, Δt = 0.001)
                # m3h = Y[:,2].^3 .* Y[:, 3] # m^3h
                # isi = collect(transpose(ISI.SpikeInterval(m3h, t, thresholdm3h)))
                # bdm3h[k] =  hcat(Astim, isi)
            end
            writedlm("results/raw/bd/bd_Amin_$(Amin)_Amax_$(Amax)_I0_$(I0)_fStim_$(fStim).csv", bd)
            #writedlm("results/raw/isi/isi_m3h_Amin_$(Amin)_Amax_$(Amax)_fStim_$(Int64(1000*fStim))_I0_$(Int64(I0)).csv", bdm3h)
        end
    end
end

#bifurcationdiag()

include("../../fhh/julia/NeuronModels.jl")

import Main.TraubMiles
import Main.Terman

function bifurcationOtherModels()
        #Constants   
    
        parameters = [#(HodgkinHuxley.SingleCell!, "HodgkinHuxley", "forced_HH", 10, 7000, 900, 0.4, 3.5, 900, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
        # (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4), 475),
                      #(TraubMiles.SingleCell!, "TraubMiles", "forced_TraubMiles", 10, 3000, 300, 0.1,10, 300, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                       #(alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                      #(WangBuzsaki.SingleCell!, "WangBuzsaki", "forced_WangBuzsaki", 10, 1000, 300, 0.1,10, 300, 0, [-6.49963792e+01, 5.95994171e-01, 3.17732402e-01],
                      # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                   #  (ErisirN2.SingleCell!, "ErisirN2", "forced_Erisir", 10, 1400, 300, 0.1,5, 300, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                    #  (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 200),
                    #  (Terman.SingleCellSTN!, "TermanSTN", "forced_TermanSTN", 10,500,200,0.1,3,200,0, [-65, 0.19, 0.15, 0.23, 0.06], 
                    # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 500),
                    #(Terman.SingleCellGP!, "TermanGP", "forced_TermanGP", 10,700,200,0.1,3,200,0, [-65, 0.19, 0.15, 0.23, 0.06], 
                    # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 500),
                    #  (RGC.typeI!, "rgcTypeI", "forced_rgcTypeI", 1, 10, 1, [1.0], [0], [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001], 
                    #  (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 80)
                    #  ]
                      # (RGC.typeII!, "rgcTypeII", "forced_rgcTypeII", 10, 600, 200, 0.2, 4.0, 200, 0, [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001], 
                      # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 100),
                      # ]
                (ReducedH3DHodgkinHuxley, "ReducedH3Dhh", "", 0, 700, 3000, 1.0, 0, [-6.483e01, 5.29550879e-02, 3.17732402e-01], (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4), 100),
                (ReducedM3DHodgkinHuxley, "ReducedM3Dhh", "", 0, 700, 3000, 1.0, 0, [-6.483e01, 5.95994171e-01, 3.17732402e-01], (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4), 100),
                (Reduced2DHodgkinHuxley, "Reduced2Dhh", "", 0, 700, 3000, 1.0, 0, [-6.483e01, 3.17732402e-01], (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4), 100)
        ]
      
         # fhh = pyimport("fhh")
      
          for (model, ModelName, pythonName, Amin, Amax, lengthA, fStimValues, I0Values, u0, diffeq, Alim) in parameters
            mkpath("results/raw/bd-$(ModelName)")
            T = 3000
            for I0 in I0Values
                for fStim in fStimValues
                    A_stimRange = collect(Amin:(Amax-Amin)/lengthA:Amax)
                    size = round(Int64, T*fStim)
                    bd = hcat(A_stimRange, zeros(length(A_stimRange), size+1))
                    l = hcat(A_stimRange, zeros(length(A_stimRange)))
                    bdm3h = Array{Union{Nothing,Matrix{Float64}}}(nothing,length(A_stimRange))
                    Threads.@threads for k in eachindex(A_stimRange)
                        Astim = A_stimRange[k]
                        println(Astim, ", ", fStim, ", ", I0)
                        p0 = [fStim, Astim, I0]
                        forcedHH = ContinuousDynamicalSystem(model.SingleCell!, u0, p0; diffeq = diffeq)
                        traj, t = trajectory(forcedHH, T; Δt = 1/fStim, Ttr = 2000)
                        bd[k, 2:end] = transpose(Matrix(traj)[:,1])
                        l[k,2] = lyapunov(forcedHH, 1000; Ttr = 500, Δt = 1)
                        if l[k,2] == NaN
                            sys = ContinuousDynamicalSystem(model.SingleCell!, u0, p0; diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4))
                            l[k,2] = lyapunov(sys, 2500; Ttr = 250, Δt = 1)[1]
                        end
                        # Y, t = trajectory(forcedHH, T; Ttr = 1000, Δt = 0.001)
                        # m3h = Y[:,2].^3 .* Y[:, 3] # m^3h
                        # isi = collect(transpose(ISI.SpikeInterval(m3h, t, thresholdm3h)))
                        # bdm3h[k] =  hcat(Astim, isi)
                    end
                    writedlm("results/raw/bd-$(ModelName)/bd_Amin_$(Amin)_Amax_$(Amax)_I0_$(I0)_fStim_$(fStim).csv", bd)
                    writedlm("results/raw/bd-$(ModelName)/lyap_Amin_$(Amin)_Amax_$(Amax)_I0_$(I0)_fStim_$(fStim).csv", l)
                    #writedlm("results/raw/isi-rgcTypeI/isi_m3h_Amin_$(Amin)_Amax_$(Amax)_fStim_$(Int64(1000*fStim))_I0_$(Int64(I0)).csv", bdm3h)
                end
            end
        end
end

bifurcationOtherModels()