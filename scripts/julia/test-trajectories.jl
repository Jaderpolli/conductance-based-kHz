include("../../fhh/julia/forcedHH.jl")
import Main.HodgkinHuxley
using DynamicalSystemsBase, Plots, OrdinaryDiffEq, LaTeXStrings, ChaosTools

function hh()

    fstim = 2.75
    Astim = 6300
    I0 = 0

    Cm  =   1.0 # membrane capacitance, in uF/cm^2
    g_Na = 120.0 # maximum conducances, in mS/cm^2
    g_K  =  36.0
    g_L  =   0.3
    E_Na =  50.0 # Nernst reversal potentials, in mV
    E_K  = -77.0
    E_L  = -54.387


    p = (fstim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L)

    u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
    diffeq = (abstol = 1e-9, reltol = 1e-9)
    forcedHH = CoupledODEs(HodgkinHuxley.SingleCell!, u0, p; diffeq)

    T = 10/fstim #ms
    Y, ty = trajectory(forcedHH, T; Ttr = 1000, Δt = 1/fstim)
    #lyap = lyapunovspectrum(forcedHH, 500)
    #println(lyap)
    #plt = plot(ty, Y[:,2].^3 .* Y[:,3], label = false, lc = :blue)
    plt = scatter(ty, Y[:, 1], lc = :red, label = false)
    #plt = plot!([1000;1030], [0.12;0.12], ls = :dot, label = false)
    # X, tx = trajectory(forcedHH, T; Ttr = 1000, Δt = 0.01)
    # plt = plot!(twinx(), [(ty, Y[:,1]), (tx, X[:,1])], ls = [:solid,:dash], lc = [:red, :yellow])
    display(plt)
end

#hh()

include("../../fhh/julia/NeuronModels.jl")
import Main.TraubMiles
import Main.ReducedTraubMiles
import Main.Erisir
import Main.WangBuzsaki
import Main.Terman
import Main.RGC

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

function trajtests()
    
    I0 = 0.0
    Astim = 6600
    fstim = 4.0 
    u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01]
    diffeq =  (alg = Tsit5(), abstol = 1e-7, reltol = 1e-7)
    #Astim = 2*pi*fstim*A
    println(Astim/(2*pi*fstim))

    #p0 = (fstim, 0, I0)
    p = (fstim, Astim, I0)

    #u0 = [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001]
    #u0 = [-5.3e01, 0.037762, 0.1, 0.018825] #4D 
    #u0=[-65, 0.19, 0.15, 0.23, 0.06] # GPe
    #u0 = [-65, 0.19, 0.15, 0.23, 0.06] # STN
   # u0 = [-65, 0.19, 0.15] #3D
    forcedHH = ContinuousDynamicalSystem(HodgkinHuxley.SingleCell!, u0, p; diffeq = diffeq)
    Y, ty = trajectory(forcedHH, 100; Ttr = 1500, Δt = 100)
    #u0 = Y[end]
    #p = (fstim, Astim, I0)
    #forcedHH = ContinuousDynamicalSystem(WangBuzsaki.SingleCell!, u0, p; diffeq = diffeq)
    lyaps = lyapunov(forcedHH, 10; Ttr = 0, Δt = 0.5)
    println(lyaps[1])
    ##maxP = 200
    #println(u0)
    #Y, ty = trajectory(forcedHH, 100; Ttr = 0, Δt = 0.01)

    #p1 = plot(ty, Y[:,1], lc = :red, label = false, frame = :box, lw = 0.5)

    # u0 = [-65, 0.59, 0.31] #3D
    # forcedHH = ContinuousDynamicalSystem(WangBuzsaki.SingleCell!, u0, p0; diffeq = diffeq)
    # Y, ty = trajectory(forcedHH, 100; Ttr = 1500, Δt = 100)
    # u0 = Y[end]
    # p = (fstim, Astim, I0)
    # forcedHH = ContinuousDynamicalSystem(WangBuzsaki.SingleCell!, u0, p; diffeq = diffeq)
    # lyaps = lyapunovspectrum(forcedHH, 100; Ttr = 0, Δt = 1)
    # println(lyaps[1])
    # maxP = 200
    # println(u0)
    # Y, ty = trajectory(forcedHH, 100; Ttr = 0, Δt = 0.01)

    #p2 = plot(ty, Y[:,1], lc = :red, label = false, frame = :box, lw = 0.5)

   # p = plot(p1, p2, layout = (2,1), plot_title = L"f_{\mathrm{stim}} = "*string(fstim)*"kHz")
    #display(p1)

    #T = 1000 #ms
    #Y, ty = trajectory(forcedHH, T; Ttr = 0, Δt = 0.01)
    #u0 = Y[end]
    #p = (fstim, Astim, I0)

    #u0 = [-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001]
    #u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01] #4D
    #u0=[-65, 0.2, 0.78, 0.9, 0.035] # GPe
    #u0 = [-65, 0.19, 0.15, 0.23, 0.06] # STN
    #u0 = [-6.65e01, 0.98633, 0.078825] #3D
    #T = 2000

    #diffeq = (alg = Rosenbrock23(), abstol = 1e-8, reltol = 1e-4)
   # forcedHH = ContinuousDynamicalSystem(RGC.typeI!, u0, p; diffeq = diffeq)
    

    #P = period(Y[:,1], maxP, 0.1)
    #println(P)
    #X = vcat(Y[:,1],Z[:,1])
    #mh = vcat(Y[:,2].^3 .* Y[:,3], Z[:,2].^3 .* Z[:,3])
    #t = vcat(ty, ty[end] .+ tz)
    #println(Y[1,1] - Y[2,1])
   #lyap = lyapunov(forcedHH, 2000; Ttr = 2000, Δt = 1)[1]
   #lyaps = lyapunovspectrum(forcedHH, 2000; Ttr = 2000, Δt = 1)
    #println(lyaps)  
    #println(lyap[1])
    #plt = plot(ty, Y[:,2].^3 .* Y[:,3], label = false, lc = :blue)
    
    #plot!(tz, Astim/(2*pi*fstim) .* sin.(2*pi*fstim .* tz .+ 0.5))
    #plt2 = plot(t, mh)
    #plt = plot!([1000;1030], [0.12;0.12], ls = :dot, label = false)
    # X, tx = trajectory(forcedHH, T; Ttr = 1000, Δt = 0.01)
    # plt = plot!(twinx(), [(ty, Y[:,1]), (tx, X[:,1])], ls = [:solid,:dash], lc = [:red, :yellow])
    
end

#trajtests()

import Main.StochasticHH
using StochasticDiffEq, DiffEqNoiseProcess, DynamicalSystems, ChaosTools

function SDEtrajTests()
    fstim = 3
    Astim = 4000
    I0 = 0

    p0 = (fstim, Astim, I0)
    u0 = [-6.49963792e+01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01] #4d
    
    
    diffeq = (alg = ImplicitEM(), abstol = 1e-3, reltol = 1e-3)
    forcedHH = CoupledSDEs(StochasticHH.SingleCell!, u0, p0; diffeq = diffeq, noise_strength = 1, covariance = [1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]);

    Y0, ty0 = trajectory(forcedHH, 10000; Ttr = 1000, Δt = 0.01)
    

    V = Y0[:,1]
    m = Y0[:,2]
    h = Y0[:,3]
    t = ty0

    lyap = lyapunov_from_data(Y0, 1000; ntype = NeighborNumber(1), distance = Euclidean())
    println(lyap[1])

    #p1 = plot(t[1:10000], m[1:10000].^3 .* h[1:10000], lc = :red, label = false, frame = :box, lw = 0.5)
    p1 = plot(t,V, lc = :red, label = false, frame = :box, lw = 0.5)
    #scatter!(t[1:100],V[101:201])
    #p2 = plot(t, (1 ./ (1 .+ exp.( .-(V .+ 30) ./ 15.0))).^3 .* h , lc = :blue, label = false, frame = :box, lw = 0.5, dpi = 300)
    #p2 = plot(t, m.^3 .* h, lc = :blue, label = false, frame = :box, lw = 0.5, dpi = 300)
    #p = plot(p1,p2,layout = (2,1))
    display(p1)   
end

#SDEtrajTests()