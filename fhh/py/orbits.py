import numpy as np
from .forcedHH import run_forced_HH, run_reduced_forced_HH, run_reduced_fast_forced_HH, run_reduced_3D_forced_HH
from .TraubMiles import run_forced_TraubMiles
from .WangBuzsaki import run_forced_WangBuzsaki
from .Erisir import run_forced_Erisir
from .TermanSTN import run_forced_TermanSTN
from .TermanGPe import run_forced_TermanGPe
from .Rgc1 import run_forced_RGC1
from .Rgc2 import run_forced_RGC2
from scipy import signal
from scipy.optimize import curve_fit

def get_orbit(f_stim, Amp, t_sim=1000, t_offset=1000, I_0=0, reduced=False, which="2D"):
    #dt = f_stim / 1e3
    dt = 1 / (f_stim * 1e3)
    if reduced is True:
        if which == "fast":
            t, X = run_reduced_fast_forced_HH(
                t_sim + t_offset,
                dt,
                I_0,
                f_stim,
                Amp,
            )
        elif which == "3D":
            t, X = run_reduced_3D_forced_HH(
                t_sim + t_offset,
                dt,
                I_0,
                f_stim,
                Amp,
            )
        else:
            t, X = run_reduced_forced_HH(
                t_sim + t_offset,
                dt,
                I_0,
                f_stim,
                Amp,
            )
    else:
        t, X = run_forced_HH(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    # get rid of transient
    t = t[(int(t_offset / dt)) :]
    t -= t[0]
    X = X[int(t_offset / dt) :, :]
    # stromboscopic postproc
    samples = np.arange(0, int(len(t)), int(1e3))
    return X[samples, 0]

def get_orbit_othermodels(f_stim, Amp, t_sim=1000, t_offset=1000, I_0=0, model="TraubMiles"):
    dt = f_stim / 1e3
    if model == "WangBuzsaki":
        t, X = run_forced_WangBuzsaki(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    if model == "Erisir":
        t, X = run_forced_Erisir(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    if model == "TermanSTN":
        t, X = run_forced_TermanSTN(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    if model == "TermanGPe":
        t, X = run_forced_TermanGPe(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc1":
        t, X = run_forced_RGC1(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc2":
        t, X = run_forced_RGC2(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    else:
        t, X = run_forced_TraubMiles(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    # get rid of transient
    t = t[(int(t_offset / dt)) :]
    t -= t[0]
    X = X[int(t_offset / dt) :, :]
    # stromboscopic postproc
    samples = np.arange(0, int(len(t)), int(1e3/f_stim))
    return X[samples, 0]

def get_average(f_stim, Amp, t_sim=1000, t_offset=1000, I_0=0, model="TraubMiles"):
    dt = f_stim / 1e3
    if model == "HH":
        t, X = run_forced_HH(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "HH3D":
        t, X = run_reduced_3D_forced_HH(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "TermanSTN":
        t, X = run_forced_TermanSTN(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "WangBuzsaki":
        t, X = run_forced_WangBuzsaki(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Erisir":
        t, X = run_forced_Erisir(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "TermanGPe":
        t, X = run_forced_TermanGPe(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc1":
        t, X = run_forced_RGC1(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc2":
        t, X = run_forced_RGC2(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    else:
        t, X = run_forced_TraubMiles(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    # get rid of transient
    t = t[(int(t_offset / dt)) :]
    t -= t[0]
    V = X[int(t_offset / dt) :, 0]
    # filter 
    fs = 1 / dt
    Q=10
    b_notch, a_notch = signal.iirnotch(f_stim, Q, fs)
    V_filtered = signal.lfilter(
        b_notch, a_notch, V
    )
    # compute and return average
    return np.mean(V_filtered)

def get_average_fit(f_stim, Amp, t_sim=1000, t_offset=1000, I_0=0, model="TraubMiles"):
    dt = f_stim / 1e3
    if model == "HH":
        t, X = run_forced_HH(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "HH3D":
        t, X = run_reduced_3D_forced_HH(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "TermanSTN":
        t, X = run_forced_TermanSTN(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "WangBuzsaki":
        t, X = run_forced_WangBuzsaki(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Erisir":
        t, X = run_forced_Erisir(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "TermanGPe":
        t, X = run_forced_TermanGPe(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc1":
        t, X = run_forced_RGC1(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    elif model == "Rgc2":
        t, X = run_forced_RGC2(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    else:
        t, X = run_forced_TraubMiles(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    # get rid of transient
    t = t[(int(t_offset / dt)) :]
    t -= t[0]
    V = X[int(t_offset / dt) :, 0]
    # filter 

    def sine_function(x, C, D):
        return Amp/(2*np.pi*f_stim) * np.sin(2*np.pi * f_stim * x + C)+D

    V_filtered = np.zeros(len(V))
    initial_guess = [3 * np.pi/2, -70]

    # Perform the curve fitting
    params, covariance = curve_fit(sine_function, t, V, p0=initial_guess)
    C, D = params
    sine = sine_function(t, C, D)
    for k, ti in enumerate(t):
        V_filtered[k] = V[k] - sine[k]+D
    return np.mean(V_filtered)