from numba import jit
import scipy.integrate as solvers
import numpy as np


# Constants
C_m = 1.0  # membrane capacitance, in uF/cm^2
g_Na = 120.0  # maximum conducances, in mS/cm^2
g_K = 36.0
g_L = 0.3
E_Na = 50.0  # Nernst reversal potentials, in mV
E_K = -77.0
E_L = -54.387
a = -1.073
b = 0.877
n0 = 0.32
h0 = 0.45

# Rates
@jit(nopython=True)
def alpha_m(V):
    return 0.1 * (V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))


@jit(nopython=True)
def beta_m(V):
    return 4.0 * np.exp(-(V + 65.0) / 18.0)


@jit(nopython=True)
def alpha_h(V):
    return 0.07 * np.exp(-(V + 65.0) / 20.0)


@jit(nopython=True)
def beta_h(V):
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))


@jit(nopython=True)
def alpha_n(V):
    return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))


@jit(nopython=True)
def beta_n(V):
    return 0.125 * np.exp(-(V + 65) / 80.0)


# conductances
@jit(nopython=True)
def gNa(m, h):
    return g_Na * m**3 * h


@jit(nopython=True)
def gK(n):
    return g_K * n**4


#  Channel currents
@jit(nopython=True)
def I_Na(V, m, h):
    return gNa(m, h) * (V - E_Na)


@jit(nopython=True) 
def I_K(V, n):
    return gK(n) * (V - E_K)


@jit(nopython=True)
def I_L(V):
    return g_L * (V - E_L)


@jit(nopython=True)
def m_inf(V):
    return alpha_m(V) / (alpha_m(V) + beta_m(V))


@jit(nopython=True)
def h_inf(V):
    return alpha_h(V) / (alpha_h(V) + beta_h(V))


@jit(nopython=True)
def n_inf(V):
    return alpha_n(V) / (alpha_n(V) + beta_n(V))


@jit(nopython=True)
def tau_m(V):
    return 1. / (alpha_m(V) + beta_m(V))


@jit(nopython=True)
def tau_h(V):
    return 1. / (alpha_h(V) + beta_h(V))


@jit(nopython=True)
def tau_n(V):
    return 1. / (alpha_n(V) + beta_n(V))


#############################
## SIMULATION OF FORCED HH ##
#############################
def run_forced_HH(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.49963792e01, 5.29550879e-02, 5.95994171e-01, 3.17732402e-01],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, m, h, n = X

        dVdt = ((I_stim(t)) - I_Na(V, m, h) - I_K(V, n) - I_L(V)) / C_m
        dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
        dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h
        dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

        return dVdt, dmdt, dhdt, dndt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X


def run_reduced_forced_HH(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.483e01, 3.17732402e-01],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, n = X

        dVdt = (
            I_stim(t)
            - g_Na
            * m_inf(V) ** 3
            * (a * n + b)
            * (V - E_Na)
            - g_K * n**4 * (V - E_K)
            - g_L * (V - E_L)
        ) / C_m
        dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

        return dVdt, dndt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X


def run_reduced_3D_forced_HH(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.483e01, 5.29550879e-02, 3.17732402e-01],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, m, n = X

        dVdt = (
            I_stim(t)
            - g_Na
            * m ** 3
            * (a * n + b)
            * (V - E_Na)
            - g_K * n**4 * (V - E_K)
            - g_L * (V - E_L)
        ) / C_m
        dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
        dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

        return dVdt, dmdt, dndt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X

def run_reduced_3Dm_forced_HH(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.483e01, 5.95994171e-01, 3.17732402e-01],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, h, n = X

        dVdt = (
            I_stim(t)
            - g_Na
            * m_inf(V) ** 3
            * h
            * (V - E_Na)
            - g_K * n**4 * (V - E_K)
            - g_L * (V - E_L)
        ) / C_m
        dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h
        dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

        return dVdt, dhdt, dndt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X


def run_reduced_fast_forced_HH(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    X_0=[-6.6e01, 1.0e-02],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, m, = X

        dVdt = ((I_stim(t)) - g_Na * m**3 * h0 * (V - E_Na) - g_K * n0**4 * (V - E_K) - I_L(V)) / C_m
        dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m

        return dVdt, dmdt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X