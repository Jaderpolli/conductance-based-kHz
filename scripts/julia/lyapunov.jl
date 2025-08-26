include("../../fhh/julia/forcedHH.jl")

import Main.LyapunovMeasures

using DelimitedFiles
using SciMLBase: ODEProblem
using DynamicalSystems, DifferentialEquations, StaticArrays, DelimitedFiles

function rangeLyapunov()
    mkpath("results/raw/lyapunov")
    # stimulation parameters lists (element j of each lists are the limits of the j-th range)
    AminList = [200]
    AmaxList = [7000]
    fMinList = [0.4]
    fMaxList = [3.5]


    lengthA = 600
    lengthf = 600

    # testing if all lists have the same length
    if length(AminList) == length(AmaxList) == length(fMinList) == length(fMaxList)
        Nlists = length(AminList)
    else
        return("invalid size of one of the parameter min/max lists")
    end

    I0List = [5, 10]

    for I0 in I0List
        for j in 1:Nlists
            println(j)
            Amin = AminList[j]
            Amax = AmaxList[j]
            fMin = fMinList[j]
            fMax = fMaxList[j]
            lyaps = LyapunovMeasures.rangeAfLyapunov(Amin, Amax, fMin, fMax, I0, lengthA, lengthf)
            writedlm("results/raw/lyapunov/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
        end
    end
end

#rangeLyapunov()

function fixfLyapunov()
    AminList = [0, 570, 620, 1000]
    AmaxList = [700, 650, 621, 2000]
    fStim = 1 #kHz

    # testing if all lists have the same length
    if length(AminList) == length(AmaxList)
        Nlists = length(AminList)
    else
        return("invalid size of one of the parameter min/max lists")
    end

    I0List = [0.0, 5.0, 10.0]

    lengthA = 1000

    for I0 in I0List
        for j in 1:Nlists
            println(j)
            AMin = AminList[j]
            AMax = AmaxList[j]
            lyaps = LyapunovMeasures.rangeALyapunov(AMin,AMax, fStim, I0, lengthA)
            writedlm("results/raw/lyapunov/lyapunov_Amin_$(AMin)_Amax_$(AMax)_fStim_$(Int64(1000*fStim))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
        end
    end
end

############### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ############
############### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ############
### !!!!!! BEFORE RUNNING rangeLyapunov, check the number of points in each range!!! !!!! ###
############### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ############
############### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ############

#rangeLyapunov()
#fixfLyapunov()


############## other models codes ######################

function rangeAfLyapunov(prob::Function, Amin, Amax, fMin, fMax, I0, lengthA, lengthf, modName, Alim)
    # stimulation parameters
    A_stimRange = collect(Amin:round(Int64,(Amax-Amin)/lengthA):Amax)
    fStimRange = collect(fMin:(fMax-fMin)/lengthf:fMax)

    # Defining the Dynamical System
    lyaps = zeros(length(fStimRange)+1, length(A_stimRange)+1)
    
    lyaps[1, 2:end] = transpose(A_stimRange)
    lyaps[2:end, 1] = fStimRange

    # i = 1
    #u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
    u0 = [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001]
    d0 = 1e-9
    
    Cm = 1.0

    Threads.@threads for fStim in fStimRange
        println(modName, ", ", fStim)
        Threads.@threads for Astim in A_stimRange
            A = Astim/(2 * pi * fStim * Cm)
            #println(Astim)
            if A < Alim
                p0 = [fStim, Astim, I0]
                forcedHH = ContinuousDynamicalSystem(prob, u0, p0)

                # Simulating

                lyap = lyapunov(forcedHH, 500; Ttr = 500, Δt = 0.1)
                lyaps[findfirst(x -> x==fStim, fStimRange)+1,findfirst(x -> x==Astim, A_stimRange) + 1] = lyap[1]
            else
                nothing
            end
        end
        writedlm("results/raw/lyapunov-$(modName)/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
    end
    return(lyaps)
end

include("../../fhh/julia/NeuronModels.jl")
import Main.TraubMiles
import Main.ReducedTraubMiles
import Main.WangBuzsaki
import Main.ErisirN2
import Main.ReducedErisirN2
import Main.ErisirN4
import Main.ReducedErisirN4

function rangeLyapunovNewModels()

    modellist = [RGC.typeI!
    #              TraubMiles.SingleCell!,
    #             ReducedTraubMiles.SingleCell!,
    #             WangBuzsaki.SingleCell!, 
    #             ErisirN2.SingleCell!,
    #             ReducedErisirN2.SingleCell!, 
    #             ErisirN4.SingleCell!,
    #             ReducedErisirN4.SingleCell!
                ]

    modelName = ["rgcTypeI"
                #"TraubMiles", 
                # "ReducedTraubMiles", 
                # "WangBuzsaki", 
                # "ErisirN2",
                # "ReducedErisirN2", 
                # "ErisirN4", 
                # "ReducedErisirN4"
                ]

    diffeqList = [(alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-2),
                (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-3),
                ]
    Alims = [100

    i = 0

    for prob in modellist
        i += 1
        println(modelName[i])
        Alim = Alims[i]
        #diffeq = diffeqList[i]
        mkpath("results/raw/lyapunov-$(modelName[i])")
        # stimulation parameters lists (element j of each lists are the limits of the j-th range)
        AminList = [10]
        AmaxList = [1000]
        fMinList = [0.1]
        fMaxList = [1.0]


        lengthA = 20
        lengthf = 20

        # testing if all lists have the same length
        if length(AminList) == length(AmaxList) == length(fMinList) == length(fMaxList)
            Nlists = length(AminList)
        else
            return("invalid size of one of the parameter min/max lists")
        end

        I0List = [0]

        for I0 in I0List
            for j in 1:Nlists
                println(j)
                Amin = AminList[j]
                Amax = AmaxList[j]
                fMin = fMinList[j]
                fMax = fMaxList[j]
                lyaps = rangeAfLyapunov(prob, Amin, Amax, fMin, fMax, I0, lengthA, lengthf, modelName[i], Alim)
                writedlm("results/raw/lyapunov-$(modelName[i])/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
            end
        end
    end
end

rangeLyapunovNewModels()