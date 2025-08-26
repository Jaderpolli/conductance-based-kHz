module HodgkinHuxley()
    using DynamicalSystems, DifferentialEquations, StaticArrays
    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[3]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[3]
        du[4] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[4]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[4]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[4] .+ p[2] * sin.(2*pi*p[1]*t) - p[8]*u[2]^(3)*u[3]*(u[1] - p[5]) - p[9]*u[4]^4*(u[1] - p[6]) - p[10]*(u[1] - p[7]))/p[3]
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[3]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[3]
        du[4] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[4]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[4]
        return SA[du[1],du[2],du[3],du[4]]
    end
end

module ReducedM3DHodgkinHuxley
    using DynamicalSystems, DifferentialEquations, StaticArrays
    function SingleCell!(du, u, p, t)
        minf = (0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)))/(0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)) + 4.0 * exp(-(u[1] + 65.0) / 18.0))
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*minf^3*u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[2]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[2]
        du[3] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[3]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[3]
        return nothing
    end
end

module ReducedH3DHodgkinHuxley
    using DynamicalSystems, DifferentialEquations, StaticArrays
    function SingleCell!(du, u, p, t)
        minf = (0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)))/(0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)) + 4.0 * exp(-(u[1] + 65.0) / 18.0))
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^3*((a * u[3] + b))*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[3]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[3]
        return nothing
    end
end


module Reduced2DHodgkinHuxley
    using DynamicalSystems, DifferentialEquations, StaticArrays
    function SingleCell!(du, u, p, t)
        minf = (0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)))/(0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0)) + 4.0 * exp(-(u[1] + 65.0) / 18.0))
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        a = -1.073
        b = 0.877
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*minf^3*(a * u[2] + b)*(u[1] - E_Na) - g_K*u[2]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[2]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[2]
        return nothing
    end
end

module ReducedTraubMiles()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 100.0 # maximum conducances, in mS/cm^2
        g_K  =  80.0
        g_L  =   0.1
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -100.0
        E_L  = -67.0
        #fStim, Astim, I0 = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(0.32 * (u[1] + 54)/(1.0 - exp(-(u[1]+54)/4))/(0.32 * (u[1] + 54)/(1 - exp(-(u[1]+54)/4))+0.28 * (u[1]+27)/(exp((u[1]+27)/5)-1)))^3 *u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        #du[2] = 0.32 * (u[1] + 54)/(1 - exp(-(u[1]+54)/4))*(1 - u[2]) - 0.28 * (u[1]+27)/(exp((u[1]+27)/5)-1) * u[2]
        du[2] = 0.128 * exp(-(u[1]+50)/18) * (1 - u[2]) - 4/(1 + exp(-(u[1]+27)/5)) * u[2]
        du[3] = 0.032 * (u[1] + 52)/(1 - exp(-(u[1]+52)/5)) * (1 - u[3]) - 0.5 * exp(-(u[1] + 57)/40) * u[3]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 100.0 # maximum conducances, in mS/cm^2
        g_K  =  80.0
        g_L  =   0.1
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -100.0
        E_L  = -67.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(0.32 * (u[1] + 54)/(1.0 - exp(-(u[1]+54)/4))/(0.32 * (u[1] + 54)/(1 - exp(-(u[1]+54)/4))+0.28 * (u[1]+27)/(exp((u[1]+27)/5)-1)))^3 *u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        #du[2] = 0.32 * (u[1] + 54)/(1 - exp(-(u[1]+54)/4))*(1 - u[2]) - 0.28 * (u[1]+27)/(exp((u[1]+27)/5)-1) * u[2]
        du[2] = 0.128 * exp(-(u[1]+50)/18) * (1 - u[2]) - 4/(1 + exp(-(u[1]+27)/5)) * u[2]
        du[3] = 0.032 * (u[1] + 52)/(1 - exp(-(u[1]+52)/5)) * (1 - u[3]) - 0.5 * exp(-(u[1] + 57)/40) * u[3]
        return SA[du[1],du[2],du[3]]
    end
end

module TraubMiles()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 100.0 # maximum conducances, in mS/cm^2
        g_K  =  80.0
        g_L  =   0.1
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -100.0
        E_L  = -67.0
        #fStim, Astim, I0 = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.32 * (u[1] + 54.0)/(1.0 - exp(-(u[1]+54.0)/4.0))*(1.0 - u[2]) - 0.28 * (u[1]+27.0)/(exp((u[1]+27.0)/5.0)-1.0) * u[2]
        du[3] = 0.128 * exp(-(u[1]+50.0)/18.0) * (1.0 - u[3]) - 4.0/(1.0 + exp(-(u[1]+27.0)/5.0)) * u[3]
        du[4] = 0.032 * (u[1] + 52.0)/(1.0 - exp(-(u[1]+52.0)/5.0)) * (1.0 - u[4]) - 0.5 * exp(-(u[1] + 57.0)/40.0) * u[4]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 100.0 # maximum conducances, in mS/cm^2
        g_K  =  80.0
        g_L  =   0.1
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -100.0
        E_L  = -67.0
        # fStim, Astim, I0 = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.32 * (u[1] + 54)/(1 - exp(-(u[1]+54)/4))*(1 - u[2]) - 0.28 * (u[1]+27)/(exp((u[1]+27)/5)-1) * u[2]
        du[3] = 0.128 * exp(-(u[1]+50)/18) * (1 - u[3]) - 4/(1 + exp(-(u[1]+27)/5)) * u[3]
        du[4] = 0.032 * (u[1] + 52)/(1 - exp(-(u[1]+52)/5)) * (1 - u[4]) - 0.5 * exp(-(u[1] + 57)/40) * u[4]
        return SA[du[1],du[2],du[3],du[4]]
    end
end

module WangBuzsaki()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 35.0 # maximum conducances, in mS/cm^2
        g_K  =  9.0
        g_L  =   0.1
        E_Na =  55.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -65.0
        #fStim, Astim, I0 = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*((0.1 * (u[1] + 35.0)/(1.0 - exp(-(u[1]+35.0)/10.0)))/(0.1 * (u[1] + 35.0)/(1.0 - exp(-(u[1]+35.0)/10.0)) + 
        4.0 * exp(-(u[1] + 60.0)/18.0)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.35 * exp(-(u[1]+58.0)/20.0) * (1.0 - u[2]) - 5.0 / (1.0 + exp(-0.1*(u[1]+28))) * u[2]
        du[3] = 0.05 * (u[1] + 34.0)/(1.0 - exp(-(u[1]+34.0)/10.0)) * (1.0 - u[3]) - 0.625 * exp(-(u[1] + 44.0)/80.0) * u[3]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 35.0 # maximum conducances, in mS/cm^2
        g_K  =  9.0
        g_L  =   0.1
        E_Na =  55.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -65.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*((0.1 * (u[1] + 35.0)/(1.0 - exp(-(u[1]+35.0)/10.0)))/(0.1 * (u[1] + 35.0)/(1.0 - exp(-(u[1]+35.0)/10.0)) + 
        4.0 * exp(-(u[1] + 60.0)/18.0)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.35 * exp(-(u[1]+58.0)/20.0) * (1.0 - u[2]) - 5.0 / (1.0 + exp(-0.1*(u[1]+28))) * u[2]
        du[3] = 0.05 * (u[1] + 34.0)/(1.0 - exp(-(u[1]+34.0)/10.0)) * (1.0 - u[3]) - 0.625 * exp(-(u[1] + 44.0)/80.0) * u[3]
        return SA[du[1],du[2],du[3]]
    end
end


module ReducedErisirN4()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        #fStim, Astim,I0  = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))))/((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) 
        + 1.2262 * exp(-u[1]/42.248)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.0035 * exp(-u[1]/24.186) * (1 - u[2]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[2]
        du[3] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[3]) - 0.025 * exp(-u[1]/22.222) * u[3]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))))/((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) 
        + 1.2262 * exp(-u[1]/42.248)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.0035 * exp(-u[1]/24.186) * (1 - u[2]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[2]
        du[3] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1.0) * (1.0 - u[3]) - 0.025 * exp(-u[1]/22.222) * u[3]
        return SA[du[1],du[2],du[3]]
    end
end

module ErisirN4()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = (40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) * (1 - u[2]) - 1.2262 * exp(-u[1]/42.248) * u[2]
        du[3] = 0.0035 * exp(-u[1]/24.186) * (1 - u[3]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[3]
        du[4] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[4]) - 0.025 * exp(-u[1]/22.222) * u[4]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = (40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) * (1 - u[2]) - 1.2262 * exp(-u[1]/42.248) * u[2]
        du[3] = 0.0035 * exp(-u[1]/24.186) * (1 - u[3]) - 0.017*(u[1] + 51.25)/(exp((-u[1]+51.25)/5.2) - 1) * u[3]
        du[4] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[4]) - 0.025 * exp(-u[1]/22.222) * u[4]        
        return SA[du[1],du[2],du[3],du[4]]
    end
end

module ReducedErisirN2()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        #fStim, Astim,I0  = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))))/((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) 
        + 1.2262 * exp(-u[1]/42.248)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^2*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.0035 * exp(-u[1]/24.186) * (1 - u[2]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[2]
        du[3] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[3]) - 0.025 * exp(-u[1]/22.222) * u[3]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*(((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))))/((40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) 
        + 1.2262 * exp(-u[1]/42.248)))^(3)*u[2]*(u[1] - E_Na) - g_K*u[3]^2*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.0035 * exp(-u[1]/24.186) * (1 - u[2]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[2]
        du[3] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1.0) * (1.0 - u[3]) - 0.025 * exp(-u[1]/22.222) * u[3]
        return SA[du[1],du[2],du[3]]
    end
end

module ErisirN2()
    using DynamicalSystems, DifferentialEquations, StaticArrays

    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^2*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = (40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) * (1 - u[2]) - 1.2262 * exp(-u[1]/42.248) * u[2]
        du[3] = 0.0035 * exp(-u[1]/24.186) * (1 - u[3]) - (-0.017)*(u[1] + 51.25)/(exp(-(u[1]+51.25)/5.2) - 1) * u[3]
        du[4] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[4]) - 0.025 * exp(-u[1]/22.222) * u[4]
        return nothing
    end

    function SingleCellStat(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 112.0 # maximum conducances, in mS/cm^2
        g_K  =  224.0
        g_L  =   0.5
        E_Na =  60.0 # Nernst reversal potentials, in mV
        E_K  = -90.0
        E_L  = -70.0
        # fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^2*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = (40 * (-u[1] + 75.5)/(exp((-u[1]+75.5)/13.5))) * (1 - u[2]) - 1.2262 * exp(-u[1]/42.248) * u[2]
        du[3] = 0.0035 * exp(-u[1]/24.186) * (1 - u[3]) - 0.017*(u[1] + 51.25)/(exp((-u[1]+51.25)/5.2) - 1) * u[3]
        du[4] = (95 - u[1])/(exp((95 - u[1])/11.8) - 1) * (1 - u[4]) - 0.025 * exp(-u[1]/22.222) * u[4]        
        return SA[du[1],du[2],du[3],du[4]]
    end
end

module Terman()
    function minf(v,thetam,sigmam)
        return 1.0/(1.0 + exp(-(v-thetam)/sigmam))
    end

    function hinf(v,thetah,sigmah)
        return 1.0/(1.0 + exp(-(v-thetah)/sigmah))
    end

    function ninf(v,thetan, sigman)
        return 1.0/(1.0 + exp(-(v-thetan)/sigman))
    end

    function ainf(v,thetaa, sigmaa)
        return 1.0/(1.0 + exp(-(v-thetaa)/sigmaa))
    end

    function rinf(v,thetar, sigmar)
        return 1.0/(1.0 + exp(-(v-thetar)/sigmar))
    end

    function sinf(v,thetas, sigmas)
        return 1.0/(1.0 + exp(-(v-thetas)/sigmas))
    end

    function binf(r,thetab,sigmab)
        return 1.0/(1.0 + exp((r - thetab)/sigmab)) - 1.0/(1.0 + exp(-thetab/sigmab))
    end

    function taun(v, taun0, taun1, thetantau, sigmantau)
        return taun0 + taun1/(1.0 + exp(-(v - thetantau)/sigmantau))
    end

    function tauh(v, tauh0, tauh1, thetahtau, sigmahtau)
        return tauh0 + tauh1/(1.0 + exp(-(v - thetahtau)/sigmahtau))
    end

    function taur(v, taur0, taur1, thetartau, sigmartau)
        return taur0 + taur1/(1.0 + exp(-(v - thetartau)/sigmartau))
    end

    function Hin(v, thetaHg, sigmaHg)
        return 1.0/(1.0 + exp(-(v - thetaHg)/sigmaHg))
    end

    function SingleCellSTN!(du, u, p, t)
            
            C_m = 1.0  # membrane capacitance, in uF/cm^2
            g_L, g_Na, g_K, g_AHP, g_Ca, g_T = (
                2.25,
                37.5,
                45.0,
                9.0,
                0.5,
                0.5,
            )  # maximum conducances, in mS/cm^2
            E_L, E_Na, E_K, E_Ca = -60, 55.0, -80.0, 140.0
            thetam, sigmam, thetah, sigmah = -30, 15.0, -39.0, -3.1
            thetan, sigman, thetar, sigmar = -32.0, 8.0, -67.0, -2.0
            thetaa, sigmaa, thetab, sigmab = -63.0, 7.8, 0.4, -0.1
            thetas, sigmas = -39.0, 8.0
            tauh0, tauh1, thetahtau, sigmahtau, ϕn, ϕh, ϕr = 1.0, 500.0, -57.0, -3.0, 0.75, 0.75, 0.2
            taun0, taun1, thetantau, sigmantau = 1.0, 100.0, -80.0, -26.0
            taur0, taur1, thetartau, sigmartau = 40.0, 17.5, 68.0, -2.2
            k1, ϵ = 15.0, 3.75e-05
            k_Ca = 22.5

            # parameters to use in network!
            # alpha, beta, thetag, gGtoS, vGtoS = 5, 1.0, 30.0, 2.25, -85
            # thetH, sigmH = -39, -8 

            du[1] = (-g_L * (u[1] - E_L) - g_K * u[3]^4 * (u[1] - E_K) - g_Na * (1.0/(1.0 + exp(-(u[1]-thetam)/sigmam)))^3 * u[2] * (u[1] - E_Na) - g_T * (1.0/(1.0 + exp(-(u[1]-thetaa)/sigmaa)))^3 * (1.0/(1.0 + exp((u[4] - thetab)/sigmab)) - 1.0/(1.0 + exp(-thetab/sigmab)))^2 * (u[1] - E_Ca) - g_Ca * (1.0/(1.0 + exp(-(u[1]-thetas)/sigmas)))^2 * (u[1] - E_Ca) - g_AHP * (u[1] - E_K) * u[5]/(u[5] + k1) + p[3] + p[2] * sin.(2 * pi * p[1] * t))/C_m
            du[2] = ϕh * ((1.0/(1.0 + exp(-(u[1]-thetah)/sigmah))) - u[2])/(tauh0 + tauh1/(1.0 + exp(-(u[1] - thetahtau)/sigmahtau)))
            du[3] = ϕn * (1.0/(1.0 + exp(-(u[1]-thetan)/sigman)) - u[3])/(taun0 + taun1/(1.0 + exp(-(u[1] - thetantau)/sigmantau)))
            du[4] = ϕr * (1.0/(1.0 + exp(-(u[1]-thetar)/sigmar)) - u[4])/(taur0 + taur1/(1.0 + exp(-(u[1] - thetartau)/sigmartau)))
            du[5] = ϵ*(-g_Ca * (1.0/(1.0 + exp(-(u[1]-thetas)/sigmas)))^2 * (u[1] - E_Ca) - g_T * (1.0/(1.0 + exp(-(u[1]-thetaa)/sigmaa)))^3 * (1.0/(1.0 + exp((u[4] - thetab)/sigmab)) - 1.0/(1.0 + exp(-thetab/sigmab)))^2 * (u[1] - E_Ca) - k_Ca * u[5])
        return nothing
    end

    function SingleCellGP!(du, u, p, t)
            
            C_m = 1.0  # membrane capacitance, in uF/cm^2
            g_L, g_Na, g_K, g_AHP, g_Ca, g_T = (
                0.1,
                120.0,
                30.0,
                30.0,
                0.15,
                0.5,
            )  # maximum conducances, in mS/cm^2
            E_L, E_Na, E_K, E_Ca = -55.0, 55.0, -80.0, 120.0
            thetam, sigmam, thetah, sigmah = -37.0, 10.0, -58.0, -12.0
            thetan, sigman, thetar, sigmar = -50.0, 14.0, -70.0, -2.0
            thetaa, sigmaa = -57.0, 2.0
            thetas, sigmas = -35.0, 2.0
            tauh0, tauh1, thetahtau, sigmahtau, ϕn, ϕh, ϕr = 0.05, 0.27, -40.0, -12.0, 0.05, 0.05, 1.0
            taun0, taun1, thetantau, sigmantau = 0.05, 0.27, -40.0, -12.0
            taurc= 30.0
            k1, ϵ = 30.0, 1e-04
            k_Ca = 20.0

            # parameters to use in network!
            # alpha, beta, thetag, gGtoS, vGtoS = 5, 1.0, 30.0, 2.25, -85
            # thetH, sigmH = -39, -8 

            du[1] = (-g_L * (u[1] - E_L) - g_K * u[2]^4 * (u[1] - E_K) - g_Na * minf(u[1],thetam, sigmam)^3 * u[3] * (u[1] - E_Na) - g_T * ainf(u[1],thetaa, sigmaa)^3 * u[4] * (u[1] - E_Ca) - g_Ca * sinf(u[1], thetas, sigmas)^2 * (u[1] - E_Ca) - g_AHP * (u[1] - E_K) * u[5]/(u[5] + k1) + p[3] + p[2] * sin.(2 * pi * p[1] * t))/C_m
            du[2] = ϕn * (ninf(u[1], thetan, sigman) - u[2])/taun(u[1],taun0, taun1, thetantau, sigmantau)
            du[3] = ϕh * (hinf(u[1], thetah, sigmah) - u[3])/tauh(u[1],tauh0, tauh1, thetahtau, sigmahtau)
            du[4] = ϕr * (rinf(u[1], thetar, sigmar) - u[4])/taurc
            du[5] = ϵ*(-g_Ca * sinf(u[1], thetas, sigmas)^2 * (u[1] - E_Ca) - g_T * ainf(u[1],thetaa, sigmaa)^3 * u[4] * (u[1] - E_Ca) - k_Ca * u[5])
        return nothing
    end
end

module RGC

    function Heavyside(x, a)
        if x > a
            return 1
        elseif x == a
            return 1/2
        elseif x < a
            return 0
        end
    end

    am(V) = -2.725 * (V + 35)/(exp(-0.1 * (V + 35)) - 1)

    bm(V) = 90.83 * exp(-(V +  60)/20)

    ah(V) = 1.817 * exp(-(V + 52)/20)

    bh(V) = 27.25 / (1 + exp(-0.1*(V+22)))

    an(V) = -0.09575 * (V + 37)/(exp(-0.1 * (V + 37)) - 1)

    bn(V) = 1.915 * exp(-(V+47)/80)

    ac(V) = -1.362 * (V + 13) / (exp(-0.1 * (V + 13)) - 1)

    bc(V) = 45.41 * exp(-(V + 38)/18)

    function typeI!(du, u, p, t)
        g_Na = 72.0 # mS/cm^2 #typeI quantity
        E_Na = 60.6 # mV
        g_L = 0.2 # mS/cm^2
        E_L = -65.02 #mV
        g_K = 50.4 #mS/cm^2  #typeI quantity
        E_K = -101.34 # mV
        g_KCa = 0.025 # ms/cm^2 (or 0.05 in Fohlmeister 1997)
        Ca_dis = 1e-6 # M
        Ca_r = 1e-7 # M
        Ca_e = 0.0018 # M
        F = 96489 # C
        r = 0.1 # um
        tauCa = 1.5 # ms
        g_Ca = 1.2 # mS/cm^2 #typeI quantity
        C_m = 1.0 # uF/cm^2
        R = 8.31451 # J/(M K)
        T = 308.15 # K (35 Celsius)
        E_Ca = R * T / (2 * F) * log(abs.(Ca_e/u[6]))

        # Dynamics
        du[1] = (-g_Na * u[2]^3 * u[3] *(u[1]- E_Na) - g_K * u[4]^4 * (u[1] - E_K) - g_KCa * (u[6]/Ca_dis)^2/(1+(u[6]/Ca_dis)^2) * (u[1] - E_K) - g_Ca * u[5]^3 * (u[1] - E_Ca) - g_L * (u[1] - E_L) + p[3] + p[2] * sin.(2 * pi * p[1] * t))/C_m
        du[2] = am(u[1]) * (1 - u[2]) - bm(u[1]) * u[2] #m
        du[3] = ah(u[1]) * (1 - u[3]) - bh(u[1]) * u[3] #h
        du[4] = an(u[1]) * (1 - u[4]) - bn(u[1]) * u[4] #n
        du[5] = ac(u[1]) * (1 - u[5]) - bc(u[1]) * u[5] #c
        du[6] = -3/(2*F*r) * g_Ca * u[5]^3 * (u[1] - E_Ca) * Heavyside(-3/(2*F*r) * g_Ca * u[5]^3 * (u[1] - E_Ca), 0) - (u[6] - Ca_r)/tauCa #Ca
        return nothing
    end

    function typeII!(du, u, p, t)
        g_Na = 88.4 # mS/cm^2 #typeII quantity
        E_Na = 60.6 # mV
        g_L = 0.2 # mS/cm^2
        E_L = -65.02 #mV
        g_K = 94.8 #mS/cm^2  #typeII quantity
        E_K = -101.34 # mV
        g_KCa = 0.025 # ms/cm^2 (or 0.05 in Fohlmeister 1997)
        Ca_dis = 1e-6 # M
        Ca_r = 1e-7 # M
        Ca_e = 0.0018 # M
        F = 96489 # C
        r = 0.1 # um
        tauCa = 1.5 # ms
        g_Ca = 0.765 # mS/cm^2 #typeII quantity
        C_m = 1.0 # uF/cm^2
        R = 8.31451 # J/(M K)
        T = 308.15 # K (35 Celsius)
        E_Ca = R * T / (2 * F) * log(abs.(Ca_e/u[6]))

        # Dynamics
        du[1] = (-g_Na * u[2]^3 * u[3] *(u[1]- E_Na) - g_K * u[4]^4 * (u[1] - E_K) - g_KCa * (u[6]/Ca_dis)^2/(1+(u[6]/Ca_dis)^2) * (u[1] - E_K) - g_Ca * u[5]^3 * (u[1] - E_Ca) - g_L * (u[1] - E_L) + p[3] + p[2] * sin.(2 * pi * p[1] * t))/C_m
        du[2] = am(u[1]) * (1 - u[2]) - bm(u[1]) * u[2] #m
        du[3] = ah(u[1]) * (1 - u[3]) - bh(u[1]) * u[3] #h
        du[4] = an(u[1]) * (1 - u[4]) - bn(u[1]) * u[4] #n
        du[5] = ac(u[1]) * (1 - u[5]) - bc(u[1]) * u[5] #c
        du[6] = -3/(2*F*r) * g_Ca * u[5]^3 * (u[1] - E_Ca) * Heavyside(-3/(2*F*r) * g_Ca * u[5]^3 * (u[1] - E_Ca), 0) - (u[6] - Ca_r)/tauCa #Ca
        return nothing
    end
end

module StochasticHH()
    
    function SingleCell!(du, u, p, t)
        Cm  =   1.0 # membrane capacitance, in uF/cm^2
        g_Na = 120.0 # maximum conducances, in mS/cm^2
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0 # Nernst reversal potentials, in mV
        E_K  = -77.0
        E_L  = -54.387    
        #fStim, Astim, Cm, I0, E_Na, E_K, E_L, g_Na, g_K, g_L = Tuple(p)
        du[1] = (p[3] .+ p[2] * sin.(2*pi*p[1]*t) - g_Na*u[2]^(3)*u[3]*(u[1] - E_Na) - g_K*u[4]^4*(u[1] - E_K) - g_L*(u[1] - E_L))/Cm
        du[2] = 0.1 * (u[1] + 40.0)/(1.0 - exp(-(u[1] + 40.0) / 10.0))* (1.0 - u[2]) - 4.0 * exp(-(u[1] + 65.0) / 18.0) * u[2]
        du[3] = 0.07 * exp(-(u[1] + 65.0) / 20.0) * (1.0 - u[3]) - 1.0 / (1 + exp(-(u[1] + 35.0) / 10.0)) * u[3]
        du[4] = 0.01 * (u[1] + 55.0) / (1.0 - exp(-(u[1] + 55.0) / 10.0)) * (1.0 - u[4]) - 0.125 * exp(-(u[1] + 65) / 80.0) * u[4]
        return nothing
    end
end