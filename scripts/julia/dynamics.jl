include("../../fhh/julia/NeuronModels.jl")

using DynamicalSystems, DelimitedFiles, DifferentialEquations#, BenchmarkTools

import Main.TraubMiles
import Main.ErisirN2
import Main.WangBuzsaki
import Main.HodgkinHuxley
import Main.Terman

function is_close(p1, p2, tol)
    return abs(p1 - p2) < tol
end

function period(X, maxP, ϵ)
    n = length(X)  # Number of points in the trajectory
    
    for P in 1:maxP
        periodic = true
        
        # Check if the trajectory repeats itself after a period P
        for i in 1:(n - P)
            if !is_close(X[i], X[i + P], ϵ)
                periodic = false
                break
            end
        end
        
        # If all points match, we found the period
        if periodic
            return P
        end
    end
    
    return 0  # Return nothing if no period <= P_max is found
end

function main()
    #= (model, ModelName, pythonName, Amin, Amax, lengthA, fMin, fMax, lengthf, I0, u0, diffeq, Alim)
        where the pythonName can be either 
            "forced_HH"
            "reduced_forced_HH"
            "reduced_3D_forced_HH"
            "reduced_fast_forced_HH"
            "forced_TraubMiles"
            "forced_WangBuzsaki"
            "forced_Erisir"
            "forced_TermanSTN"
            "forced_TermanGPe"
    =#

    parameters = [#(HodgkinHuxley.SingleCell!, "HodgkinHuxley", "forced_HH", 10, 7000, 20, 2.5, 6.0, 20, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
   #(alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4), 425),
                (StochasticHH.SingleCell!, "Stoch_HodgkinHuxley", "forced_HH", 10, 7000, 20, 0.4, 6.0, 20, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                 (alg = ImplicitEM(), abstol = 1e-3, reltol = 1e-3), 425),
                # (TraubMiles.SingleCell!, "TraubMiles", "forced_TraubMiles", 10, 3000, 300, 0.1, 7.5, 300, 10, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                #(TraubMiles.SingleCell!, "TraubMiles", "forced_TraubMiles", 10, 4000, 100, 0.1, 7.5, 100, 0, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                #(alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                # (TraubMiles.SingleCell!, "TraubMiles", "forced_TraubMiles", 10, 4000, 100, 0.1, 7.5, 100, 20, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                # (WangBuzsaki.SingleCell!, "WangBuzsaki", "forced_WangBuzsaki", 10, 500, 500, 0.1, 5, 500, 0, [-6.49963792e+01, 5.95994171e-01, 3.17732402e-01],
               # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                  #(WangBuzsaki.SingleCell!, "WangBuzsaki", "forced_WangBuzsaki", 10, 500, 500, 0.1,5, 500, 10, [-6.49963792e+01, 5.95994171e-01, 3.17732402e-01],
                  #(alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 120),
               #   (WangBuzsaki.SingleCell!, "WangBuzsaki", "forced_WangBuzsaki", 50, 150, 1000, 0.7, 1.0, 1000, 10, [-6.49963792e+01, 5.95994171e-01, 3.17732402e-01],
               #   (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                #(ErisirN2.SingleCell!, "ErisirN2", "forced_Erisir", 10, 1400, 300, 0.1,5, 300, 10, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                #(ErisirN2.SingleCell!, "ErisirN2", "forced_Erisir", 10, 1400, 900, 0.1,5, 900, 10, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                #(alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
                #(ErisirN2.SingleCell!, "ErisirN2", "forced_Erisir", 10, 1400, 900, 0.1,5, 900, 20, [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
                #(alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 150),
               #  (Terman.SingleCellSTN!, "TermanSTN", "forced_TermanSTN", 10,1500,400,0.1,10,400,-10, [-65, 0.19, 0.15, 0.23, 0.06], 
               # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 500),
           # (Terman.SingleCellGP!, "TermanGP", "forced_TermanGP", 10,1500,300,0.8, 5, 300,-5, [-65, 0.19, 0.15, 0.23, 0.06], 
            # (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 500),
            #   (RGC.typeI!, "rgcTypeI", "forced_rgcTypeI", 10, 4000, 200, 0.5, 15, 200, -5, [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001], 
            #  (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 100),
              #(RGC.typeII!, "rgcTypeII", "forced_rgcTypeII", 10, 4000, 200, 0.5, 15, 200, -5, [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001], 
             #   (alg = Tsit5(), abstol = 1e-6, reltol = 1e-6), 100)
               ]

   # fhh = pyimport("fhh")

    for (model, ModelName, pythonName, Amin, Amax, lengthA, fMin, fMax, lengthf, I0, u0, diffeq, Alim) in parameters
        mkpath("results/raw/dynamics-$(ModelName)")
        A_stimRange = collect(Amin:(Amax-Amin)/lengthA:Amax)
        fStimRange = collect(fMin:(fMax-fMin)/lengthf:fMax)

       # creating lyapunov matrix
        lyaps = zeros(length(fStimRange)+1, length(A_stimRange)+1)
        lyaps[1, 2:end] = transpose(A_stimRange)
        lyaps[2:end, 1] = fStimRange

        # creating period matrix
        M = (-2) .* ones(length(fStimRange)+1, length(A_stimRange)+1)
        M[1, 2:end] = transpose(A_stimRange)
        M[2:end, 1] = fStimRange

        threshold = 1e-1
        maxP = 200 
        

        if ModelName == "Stoch_HodgkinHuxley"
            sys1 = CoupledSDEs(model, u0, [1,0,I0]; diffeq = diffeq, noise_strength = 1.0, covariance = [1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])  
            traj1, ttraj =  trajectory(sys1, 1000; Ttr = 100, Δt = 0.01)
            u0 = traj1[end] 
            Threads.@threads for i in eachindex(fStimRange)
                fStim = fStimRange[i]
                println(ModelName, ", ", fStim)
                for j in eachindex(A_stimRange)
                    Astim = A_stimRange[j]
                    A = Astim/(2*pi*fStim)
                    if A < Alim
                        sys2 = CoupledSDEs(model, u0, [fStim,Astim,I0]; diffeq = diffeq, noise_strength = 1.0, covariance = [1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])
                        Y0, ty0 = trajectory(sys2, 1000; Ttr = 1000, Δt = 0.01)
                        lyap = lyapunov_from_data(Y0, 100; ntype = NeighborNumber(5), distance = Euclidean())[1]
                        lyaps[i+1,j+1] = lyap
                        if lyap > threshold
                            M[i+1,j+1] = -1
                        elseif 0 < lyap < threshold
                            M[i+1,j+1] = 0
                        elseif lyap < 0
                            spike = any(Y0[:,2].^3 .* Y0[:,3] .> 0.12)
                            P = spike ? 1.0 : -3.0
                            # N = 100
                            # Y, ty = trajectory(sys2, maxP*1/fStim; Ttr = 1000, Δt = 1/(fStim))
                            # P = period(Y[:,1], maxP, 0.1)
                            # if P == 1
                            #     Y,ty = trajectory(sys2, 1000; Ttr = 1000, Δt = 0.01)
                            #     if ModelName == "TermanGP"
                            #         spike = any((1 ./ (1 .+ exp.( .-(Y[:,1]  .+ 37) ./ 10.0))).^3 .* Y[:,2] .> 0.12)
                            #     elseif ModelName == "TermanSTN"
                            #         spike = any((1 ./ (1 .+ exp.( .-(Y[:,1] .+ 30) ./ 15.0))).^3 .* Y[:,2] .> 0.12)
                            #     elseif ModelName == "WangBuzsaki"
                            #         spike = any(((0.1 .* (Y[:,1] .+ 35.0) ./ (1.0 .- exp.(.-(Y[:,1] .+ 35.0) ./ 10.0))) ./ (0.1 .* (Y[:,1] .+ 35.0)./(1.0 .- exp.(.-Y[:,1] .+ 35.0) ./ 10.0)) .+ 4.0 .* exp.(.-(Y[:,1] .+ 60.0) ./18.0)).^3 .* Y[:, 2] .> 0.12)
                            #     else
                            #         spike = any(Y[:,2].^3 .* Y[:,3] .> 0.12)
                            #     end
                            #     P = spike ? 1.0 : -3.0
                            # end                        
                            M[i+1,j+1] = P
                        end
                    end
                end
                writedlm("results/raw/dynamics-$(ModelName)/periods_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", M)
            end
            writedlm("results/raw/dynamics-$(ModelName)/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
            writedlm("results/raw/dynamics-$(ModelName)/periods_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", M)
        else
            sys1 = ContinuousDynamicalSystem(model, u0, [1, 0, I0]; diffeq = diffeq)       
            traj1, ttraj = trajectory(sys1, 1000; Ttr = 100, Δt = 0.01)
            u0 = traj1[end] 
            Threads.@threads for i in eachindex(fStimRange)
                fStim = fStimRange[i]
                println(ModelName, ", ", fStim)
                for j in eachindex(A_stimRange)
                    Astim = A_stimRange[j]
                    A = Astim/(2*pi*fStim)
                    if A < Alim
                        sys2 = ContinuousDynamicalSystem(model, u0, [fStim, Astim, I0]; diffeq = diffeq)
                            # Simulating
                        #lyap = lyapunov(sys, 2000; Ttr = 2000, Δt = 1)
                        lyap = lyapunovspectrum(sys2, 2000; Ttr = 500, Δt = 1)[1]
                        lyaps[i+1,j+1] = lyap

                        #without Simulating
                        #lyap[i,j] = lyaps[i+1,j+1]
                        # println(lyap[1])
                        if isnan(lyap)
                            sysStiff = ContinuousDynamicalSystem(model, u0,  [fStim, Astim, I0]; diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-2))
                            #lyap[1] = lyapunov(sys2, 500; Ttr = 500, Δt = 0.1)[1]
                            lyap = lyapunovspectrum(sysStiff, 500; Ttr = 500, Δt = 1)
                            lyaps[i+1,j+1] = lyap
                        end

                        if isnan(lyap)
                            M[i+1,j+1] = -2
                        elseif lyap > threshold
                            M[i+1,j+1] = -1
                        elseif 0 < lyap < threshold
                            M[i+1,j+1] = 0
                        elseif lyap < 0
                            N = 100
                            Y, ty = trajectory(sys2, maxP*1/fStim; Ttr = 1000, Δt = 1/(fStim))
                            P = period(Y[:,1], maxP, 0.1)
                            if P == 1
                                Y,ty = trajectory(sys2, 1000; Ttr = 1000, Δt = 0.01)
                                if ModelName == "TermanGP"
                                    spike = any((1 ./ (1 .+ exp.( .-(Y[:,1]  .+ 37) ./ 10.0))).^3 .* Y[:,2] .> 0.12)
                                elseif ModelName == "TermanSTN"
                                    spike = any((1 ./ (1 .+ exp.( .-(Y[:,1] .+ 30) ./ 15.0))).^3 .* Y[:,2] .> 0.12)
                                elseif ModelName == "WangBuzsaki"
                                    spike = any(((0.1 .* (Y[:,1] .+ 35.0) ./ (1.0 .- exp.(.-(Y[:,1] .+ 35.0) ./ 10.0))) ./ (0.1 .* (Y[:,1] .+ 35.0)./(1.0 .- exp.(.-Y[:,1] .+ 35.0) ./ 10.0)) .+ 4.0 .* exp.(.-(Y[:,1] .+ 60.0) ./18.0)).^3 .* Y[:, 2] .> 0.12)
                                else
                                    spike = any(Y[:,2].^3 .* Y[:,3] .> 0.12)
                                end
                                P = spike ? 1.0 : -3.0
                            end                        
                            M[i+1,j+1] = P
                        end
                    end
                end
                writedlm("dynamics-$(ModelName)/periods_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", M)
            end
            writedlm("results/raw/dynamics-$(ModelName)/lyapunov_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", lyaps)
            writedlm("results/raw/dynamics-$(ModelName)/periods_Amin_$(Amin)_Amax_$(Amax)_fMin_$(Int64(1000*fMin))_fMax_$(Int64(1000*fMax))_I0_$(Int64(I0))_size_$(lengthA)_x_$(lengthf).csv", M)
        end

        

    end
end

main()
