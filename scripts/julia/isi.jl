include("../../fhh/julia/forcedHH.jl")

import Main.ISI

using DelimitedFiles, DynamicalSystems, DifferentialEquations

function isi_bd()
    # limits
    Amin = [0, 0, 500, 1000]
    Amax = [600, 800, 2000, 3000]

    #Constants

    Cm  =   1.0 # membrane capacitance, in uF/cm^2
    g_Na = 120.0 # maximum conducances, in mS/cm^2
    g_K  =  36.0
    g_L  =   0.3
    E_Na =  50.0 # Nernst reversal potentials, in mV
    E_K  = -77.0
    E_L  = -54.387


    # stimulation parameters
    lengthA = 1000
    fStimRange = [0.5, 1.0, 1.2, 2.0]
    I0s = [0.0, 5.0, 10.0]

    #initial conditions
    u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
    # integrator
    diffeq = (abstol = 1e-8, reltol = 1e-8)
    # integration time in ms
    T = 1000 #ms
    thresholdV = 0.0
    thresholdm3h = 0.12

    mkpath("results/raw/isi")

    for I0 in I0s
        for i in eachindex(fStimRange)
            fStim = fStimRange[i]
            A_stimRange = collect(range(Amin[i], Amax[i], lengthA))
            bdV = Array{Union{Nothing,Matrix{Float64}}}(nothing,length(A_stimRange))
            bdm3h = Array{Union{Nothing,Matrix{Float64}}}(nothing,length(A_stimRange))
            Threads.@threads for j in eachindex(A_stimRange)
                Astim = A_stimRange[j]
                println(Astim)
                p0 = (fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L)
                forcedHH = CoupledODEs(HodgkinHuxley.SingleCell!, u0, p0; diffeq)
                Y, t = trajectory(forcedHH, T; Ttr = 1000, Δt = 0.001)
                V = Y[:,1]
                m3h = Y[:,2].^3 .* Y[:, 3] # m^3h
                #display(plot(t[1:10000], V[1:10000]))
                isi = collect(transpose(ISI.SpikeInterval(V, t, thresholdV)))
                #println(isi)
                bdV[j] = hcat(Astim, isi)
                isi = collect(transpose(ISI.SpikeInterval(m3h, t, thresholdm3h)))
                bdm3h[j] =  hcat(Astim, isi)
            end
            writedlm("results/raw/isi/isi_V_Amin_$(Amin[i])_Amax_$(Amax[i])_fStim_$(Int64(1000*fStim))_I0_$(Int64(I0)).csv", bdV)
            writedlm("results/raw/isi/isi_m3h_Amin_$(Amin[i])_Amax_$(Amax[i])_fStim_$(Int64(1000*fStim))_I0_$(Int64(I0)).csv", bdm3h)
        end
    end
end

isi_bd()