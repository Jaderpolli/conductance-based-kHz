from numba import jit
import scipy.integrate as solvers
import numpy as np

# Constants
C_m = 1.0  # membrane capacitance, in uF/cm^2
g_Na = 35.0  # maximum conducances, in mS/cm^2
g_K = 9.0
g_L = 0.1
E_Na = 55.0  # Nernst reversal potentials, in mV
E_K = -90.0
E_L = -65.0

# Rates
@jit(nopython=True)
def alpha_m(V):
    return 0.1 * (V + 35.0) / (1.0 - np.exp(-(V + 35.0) / 10.0))


@jit(nopython=True)
def beta_m(V):
    return 4.0 * np.exp(-(V + 60.0) / 18.0)


@jit(nopython=True)
def alpha_h(V):
    return 0.35 * np.exp(-(V + 58.0) / 20.0)


@jit(nopython=True)
def beta_h(V):
    return 5.0 / (1.0 + np.exp(-0.1 * (V + 28.0)))


@jit(nopython=True)
def alpha_n(V):
    return 0.05 * (V + 34.0) / (1.0 - np.exp(-0.1 * (V + 34.0)))


@jit(nopython=True)
def beta_n(V):
    return 0.625 * np.exp(-(V + 44.0) / 80.0)

@jit(nopython=True)
def m_inf(V):
    return alpha_m(V) / (alpha_m(V) + beta_m(V))


# conductances
@jit(nopython=True)
def gNa(h, V):  # m, 
    return g_Na * m_inf(V)**3 * h


@jit(nopython=True)
def gK(n):
    return g_K * n**4


#  Channel currents
@jit(nopython=True)
def I_Na(V, h):  # , m
    return gNa(h, V) * (V - E_Na)


@jit(nopython=True) 
def I_K(V, n):
    return gK(n) * (V - E_K)


@jit(nopython=True)
def I_L(V):
    return g_L * (V - E_L)


@jit(nopython=True)
def m_inf(V):
    return alpha_m(V) / (alpha_m(V) + beta_m(V))


#############################
## SIMULATION OF FORCED HH ##
#############################
def run_forced_WangBuzsaki(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.65e01, 0.98633, 0.078825],  # , 0.037762
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, h, n = X  # , m

        dVdt = ((I_stim(t)) - I_Na(V, h) - I_K(V, n) - I_L(V)) / C_m  # , m
        #  dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
        dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h
        dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

        return dVdt, dhdt, dndt  # , dmdt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X
