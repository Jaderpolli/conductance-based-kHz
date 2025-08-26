from numba import jit
import scipy.integrate as solvers
import numpy as np


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

#######  Rate Functions  #######
@jit(nopython=True)
def am(V):
    return -2.725 * (V + 35)/(np.exp(-0.1 * (V + 35)) - 1)

@jit(nopython=True)
def bm(V):
    return 90.83 * np.exp(-(V +  60)/20)

@jit(nopython=True)
def ah(V):
    return 1.817 * np.exp(-(V + 52)/20)

@jit(nopython=True)
def bh(V):
    return 27.25 / (1 + np.exp(-0.1*(V+22)))

@jit(nopython=True)
def an(V):
    return -0.09575 * (V + 37)/(np.exp(-0.1 * (V + 37)) - 1)

@jit(nopython=True)
def bn(V):
    return 1.915 * np.exp(-(V+47)/80)

@jit(nopython=True)
def ac(V):
    return -1.362 * (V + 13) / (np.exp(-0.1 * (V + 13)) - 1)

@jit(nopython=True)
def bc(V):
    return 45.41 * np.exp(-(V + 38)/18)

# Currents
@jit(nopython=True)
def il(v):
    return g_L * (v - E_L)


@jit(nopython=True)
def ina(v, m, h):
    return g_Na * m ** 3 * h * (v - E_Na)


@jit(nopython=True)
def ik(v, n):
    return g_K * n**4 * (v - E_K)


@jit(nopython=True)
def ikca(v, ca):
    return g_KCa * (ca/Ca_dis)**2/(1+(ca/Ca_dis)**2) * (v - E_K)


@jit(nopython=True)
def ica(v, c, E_Ca):
    return g_Ca * c**3 * (v - E_Ca)

@jit(nopython=True)
def curr(v, m, h, n, c, ca, E_Ca):
    return il(v) + ina(v, m, h) + ik(v, n) + ikca(v, ca) + ica(v, c, E_Ca)

@jit(nopython=True)
def heaviside(x1,x2):
    if x1 > x2:
        return 1.0
    elif x1 == x2:
        return 0.5
    else:
        return 0

def run_forced_RGC2(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    #  X_0=[-77, 0.19, 0.15, 0.23, 0.06, 0.0],
     X_0=[-70, 0.0353, 0.9054, 0.0677, 0.0019, 0.0000001],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, m, h, n, c, Ca = X
        E_Ca = R * T / (2 * F) * np.log(abs(Ca_e/Ca))
        H = heaviside(-3/(2*F*r) * g_Ca * c**3 * (V - E_Ca), 0.0)
        dVdt = ((I_stim(t)) - (curr(V, m, h, n, c, Ca, E_Ca))) / C_m
        dmdt = am(V) * (1 - m) - bm(V) * m #m
        dhdt = ah(V) * (1 - h) - bh(V) * h #h
        dndt = an(V) * (1 - n) - bn(V) * n #n
        dcdt = ac(V) * (1 - c) - bc(V) * c #c
        dCadt = -3/(2*F*r) * g_Ca * c**3 * (V - E_Ca) * H - (Ca - Ca_r)/tauCa #Ca

        return dVdt, dmdt, dhdt, dndt, dcdt, dCadt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X