from numba import jit
import scipy.integrate as solvers
import numpy as np

# Constants
C_m = 1.0  # membrane capacitance, in uF/cm^2
gnag, gkg, gahpg, gtg, gcag, glg = 120.0, 30.0, 30.0, 0.5, 0.15, 0.1
vnag, vkg, vcag, vlg = 55.0, -80.0, 120.0, -55.0
thmg, sigmg = -37.0, -10.0
thhg, sighg = -58.0, 12.0
thng, signg = -50.0, -14.0
thrg, sigrg, taurg = -70.0, 2.0, 30.0
thag, sigag, thsg, sigsg = -57.0, -2.0, -35.0, -2.0
tauhg0, tauhg1, thhgt, shg = 0.05, 0.27, -40, 12.0
taung0, taung1, thngt, sng = 0.05, 0.27, -40, 12.0
k1g, kcag, epsg = 30.0, 20.0, 0.0001
phirg, phing, phihg = 1.0, 0.05, 0.05
ig = -1.0
gStoG, vStoG, alphag, betag, thetagg = 0.3, 0.0, 2.0, 0.04, 20.0
vGtoG, gGtoG = -85.0, 0.1
thetHg, sigmHg = -57.0, -2.0


#######  Rate Functions  #######
@jit(nopython=True)
def minfg(v):
    return 1.0 / (1.0 + np.exp((v - thmg) / sigmg))


@jit(nopython=True)
def hinfg(v):
    return 1 / (1 + np.exp((v - thhg) / sighg))


@jit(nopython=True)
def ninfg(v):
    return 1.0 / (1.0 + np.exp((v - thng) / signg))


@jit(nopython=True)
def rinfg(v):
    return 1 / (1 + np.exp((v - thrg) / sigrg))


@jit(nopython=True)
def ainfg(v):
    return 1 / (1 + np.exp((v - thag) / sigag))


@jit(nopython=True)
def sinfg(v):
    return 1 / (1 + np.exp((v - thsg) / sigsg))


@jit(nopython=True)
def tauhg(v):
    return tauhg0 + tauhg1 / (1 + np.exp((v - thhgt) / shg))


@jit(nopython=True)
def taung(v):
    return taung0 + taung1 / (1 + np.exp((v - thngt) / sng))


####### GPe Currents  #######
@jit(nopython=True)
def ilg(v):
    return glg * (v - vlg)


@jit(nopython=True)
def inag(v, h):
    return gnag * (minfg(v) ** 3) * h * (v - vnag)


@jit(nopython=True)
def ikg(v, n):
    return gkg * (n**4) * (v - vkg)


@jit(nopython=True)
def iahpg(v, ca):
    return gahpg * (v - vkg) * ca / (ca + k1g)


@jit(nopython=True)
def icag(v):
    return gcag * ((sinfg(v)) ** 2) * (v - vcag)


@jit(nopython=True)
def itg(v, r):
    return gtg * (ainfg(v) ** 3) * r * (v - vcag)


@jit(nopython=True)
def currg(v, h, n, ca, r):
    return itg(v, r) + inag(v, h) + ikg(v, n) + ilg(v) + iahpg(v, ca) + icag(v)


def run_forced_TermanGPe(
    t_stop,
    dt,
    I_0,
    f_stim,
    Amp,
    #  X_0=[-77, 0.19, 0.15, 0.23, 0.06, 0.0],
    X_0=[-65, 0.2, 0.78, 0.9, 0.035],
):
    @jit(nopython=True)
    def I_stim(t):
        return I_0 + Amp * np.sin(2 * np.pi * f_stim * t)

    # differential equation system
    @jit(nopython=True)
    def dXdt(X, t):
        V, h, n, r, Ca = X
        dVdt = ((I_stim(t)) - (currg(V, h, n, Ca, r))) / C_m
        dndt = phing * (ninfg(V) - n) / taung(V)
        dhdt = phihg * (hinfg(V) - h) / tauhg(V)
        drdt = phirg * (rinfg(V) - r) / taurg
        dCadt = epsg * (-icag(V) - itg(V, r) - kcag * Ca)

        return dVdt, dhdt, dndt, drdt, dCadt

    t = np.linspace(0, t_stop, int(t_stop / dt))
    X, infos = solvers.odeint(dXdt, X_0, t, full_output=True, printmessg=False)
    if infos['message'] != 'Integration successful.':
        X = None

    return t, X
