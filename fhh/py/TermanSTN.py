from numba import jit
import scipy.integrate as solvers
import numpy as np

# Constants
C_m = 1.0  # membrane capacitance, in uF/cm^2
g_L, g_Na, g_K, g_AHP, g_Ca, g_T = (
    2.25,
    37.5,
    45,
    9.0,
    0.5,
    0.5,
)  # maximum conducances, in mS/cm^2
E_L, E_Na, E_K, E_Ca = -60, 55.0, -80.0, 140.0
thetam, sigmam, thetah, sigmah = -30, -15.0, -39.0, 3.1
thetan, sigman, thetar, sigmar = -32.0, -8.0, -67, 2.0
thetaa, sigmaa, thetab, sigmab = -63.0, -7.8, 0.4, -0.1
thetas, sigmas = -39, -8
tauh0, tauh1, thh, sigmaht, phi = 1, 500.0, -57.0, 3.0, 0.75
taun0, taun1, thn, sigmant = 1, 100.0, -80.0, 26.0
taur0, taur1, thr, sigmart, phir = 40, 17.5, 68, 2.2, 0.2
k1, eps = 15.0, 5e-05
k_Ca = 22.5
alpha, beta, thetag, gGtoS, vGtoS = 5, 1.0, 30.0, 2.25, -85
thetH, sigmH = -39, -8


#######  Rate Functions  #######
@jit(nopython=True)
def minf(v):
    return 1.0 / (1.0 + np.exp((v - thetam) / sigmam))


@jit(nopython=True)
def hinf(v):
    return 1.0 / (1.0 + np.exp((v - thetah) / sigmah))


@jit(nopython=True)
def ninf(v):
    return 1.0 / (1.0 + np.exp((v - thetan) / sigman))


@jit(nopython=True)
def rinf(v):
    return 1 / (1 + np.exp((v - thetar) / sigmar))


@jit(nopython=True)
def ainf(v):
    return 1 / (1 + np.exp((v - thetaa) / sigmaa))


@jit(nopython=True)
def binf(r):
    return 1 / (1 + np.exp((r - thetab) / sigmab)) - 1 / (1 + np.exp(-thetab / sigmab))


@jit(nopython=True)
def sinf(v):
    return 1.0 / (1.0 + np.exp((v - thetas) / sigmas))


@jit(nopython=True)
def tauh(v):
    return tauh0 + tauh1 / (1 + np.exp((v - thh) / sigmaht))


@jit(nopython=True)
def taun(v):
    return taun0 + taun1 / (1 + np.exp((v - thn) / sigmant))


@jit(nopython=True)
def taur(v):
    return taur0 + taur1 / (1 + np.exp((v - thr) / sigmart))


@jit(nopython=True)
def Hin(v):
    return 1 / (1 + np.exp((v - thetH) / sigmH))


####### STN Currents  #######
@jit(nopython=True)
def il(v):
    return g_L * (v - E_L)


@jit(nopython=True)
def ina(v, h):
    return g_Na * (minf(v)) ** 3 * h * (v - E_Na)


@jit(nopython=True)
def ik(v, n):
    return g_K * n**4 * (v - E_K)


@jit(nopython=True)
def iahp(v, ca):
    return g_AHP * (v - E_K) * ca / (ca + k1)


@jit(nopython=True)
def ica(v):
    return g_Ca * ((sinf(v)) ** 2) * (v - E_Ca)


@jit(nopython=True)
def it(v, r):
    return g_T * (ainf(v) ** 3) * (binf(r) ** 2) * (v - E_Ca)


@jit(nopython=True)
def curr(v, h, n, ca, r):
    return il(v) + ina(v, h) + ik(v, n) + iahp(v, ca) + ica(v) + it(v, r)


def run_forced_TermanSTN(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    #  X_0=[-77, 0.19, 0.15, 0.23, 0.06, 0.0],
    X_0=[-65, 0.19, 0.15, 0.23, 0.06],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        #V, h, n, r, Ca, s = X
        V, h, n, r, Ca = X
        dVdt = ((I_stim(t)) - (curr(V, h, n, Ca, r))) / C_m
        dhdt = phi * (hinf(V) - h) / tauh(V)
        dndt = phi * (ninf(V) - n) / taun(V)
        drdt = phir * (rinf(V) - r) / taur(V)
        dCadt = (phi * eps) * (-g_Ca * ((sinf(V))**2) * (V - E_Ca) - it(V, r) - k_Ca * Ca)
        #dsdt = alpha * (1 - s) * Hin(V - thetag) - beta * s

        return dVdt, dhdt, dndt, drdt, dCadt  # , dsdt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X
