from .forcedHH import *
from .TraubMiles import run_forced_TraubMiles
from .WangBuzsaki import run_forced_WangBuzsaki
from .Erisir import run_forced_Erisir
from .TermanSTN import run_forced_TermanSTN
from .TermanGPe import run_forced_TermanGPe

def can_you_integrate(model, freq, A, t_stop = 10, dt=0.01, I_0=0):
    """Check if the model can be integrated for the given frequency and amplitude

    Args:
        model (str): name of the model as a string from the following list:
            "forced_HH"
            "reduced_forced_HH"
            "reduced_3D_forced_HH"
            "reduced_fast_forced_HH"
            "forced_TraubMiles"
            "forced_WangBuzsaki"
            "forced_Erisir"
            "forced_TermanSTN"
            "forced_TermanGPe"
        freq (float): _description_
        A (float): _description_

    Returns:
        Flag (bool): True if the model can be integrated, False otherwise,
            if the model is not in the list, Flag is False

    Note:
        lsoda prints a warning whatsoever if the model is not integrable...
        it is going to rubish your prompt anyway
    """
    Flag = True
    if model == "forced_HH":
        t, X = run_forced_HH(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "reduced_forced_HH":
        t, X = run_reduced_forced_HH(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "reduced_3D_forced_HH":
        t, X = run_reduced_3D_forced_HH(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "reduced_fast_forced_HH":
        t, X = run_reduced_fast_forced_HH(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "forced_TraubMiles":
        t, X = run_forced_TraubMiles(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "forced_WangBuzsaki":
        t, X = run_forced_WangBuzsaki(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "forced_Erisir":
        t, X = run_forced_Erisir(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "forced_TermanSTN":
        t, X = run_forced_TermanSTN(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    elif model == "forced_TermanGPe":
        t, X = run_forced_TermanGPe(t_stop, dt, I_0, freq, A)
        if X is None:
            Flag = False
    else:
        Flag = False

    return Flag
