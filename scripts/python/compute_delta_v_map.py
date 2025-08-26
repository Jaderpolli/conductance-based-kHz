import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import sys
import fhh
from scipy import signal
from scipy.optimize import curve_fit

def sine_function(x, C, f_stim):
    return np.sin(2*np.pi * f_stim * x + C)

def filteredV1(X, dt, f_stim):
    # filter 
    fs = 1 / dt
    Q=1
    b_notch, a_notch = signal.iirnotch(f_stim, Q, fs)
    V_filtered = signal.lfilter(
        b_notch, a_notch, X
    )
    return V_filtered

def filteredV(X, t, Amp, f_stim):

    V_filtered = np.zeros(len(X))
    bounds = ([-100, f_stim, -np.pi, -np.inf], [100, f_stim, np.pi, np.inf])
    initial_guess = [3 * np.pi/2, f_stim]
    # Perform the curve fitting
    params, covariance = curve_fit(sine_function, t, X, p0=initial_guess)
    C, f_stim_fit = params
    sine = sine_function(t, C, f_stim_fit)
    for k, ti in enumerate(t):
        #sine[k] = Amp / (2*np.pi*f_stim) * np.sin(2*np.pi*f_stim * ti + 3*np.pi/2- np.pi/64)
        V_filtered[k] = X[k] - Amp/(2*np.pi*f_stim) * sine[k]
    return(V_filtered, sine)

def test_filter():
    t_offset = 1000
    I_0 = 0
    f_stim = 4.0
    Amp = 600
    dt = 0.01
    t_sim = 10
    t, X = fhh.run_forced_RGC1(
            t_sim + t_offset,
            dt,
            I_0,
            f_stim,
            Amp,
        )
    t = t[(int(t_offset / dt)) :]
    t -= t[0]
    V = X[int(t_offset / dt) :, 0]
    print(Amp/(2*np.pi*f_stim))
    #Vf = filteredV1(V, dt, f_stim)
    Vf, sine = filteredV(V, t, Amp, f_stim)
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (10,6))
    ax1.plot(t, V, c = 'black')
    ax3.plot(t, Vf, c = 'red')
    for k, s in enumerate(sine):
        sine[k] = Amp/(2*np.pi*f_stim) * sine[k] - 80
    ax2.plot(t, sine, c = 'blue')
    plt.show()

#test_filter()

# [modelname, Amin, Amax, Alen, fmin, fmax, flen, I0, Alim]
models = [
    #["WangBuzsaki", 10, 1000, 100, 0.1, 10.0, 100, 10, 200],
    #["TraubMiles", 10, 3000, 100, 0.1, 10.0, 100, 10, 200],
    #["Erisir", 10, 1400, 100, 0.1, 5.0, 100, 10, 200],
    #["TermanGPe", 10, 600, 100, 0.1, 3.0, 100, 10, 500],
    #["TermanSTN", 10, 2000, 100, 0.1, 10.0, 100, 10, 500],
    #["HH", 10, 7000, 100, 0.1, 3.5, 100, 0, 425],
    #["Rgc1", 10, 4000, 100, 0.5, 20.0, 100, 10, 100],
    ["Rgc2", 10, 4000, 100, 0.5, 20.0, 100, 10, 100]
    ]

# models = ["Rgc1", "Rgc2"]
# I_0 = 0
# fmin = 0.1
# fmax = 2.5
# nfreq = 20


Ntask = 31

for modelname, Amin, Amax, Alen, fmin, fmax, flen, I0, Alim in models:
    fname = "average_map_"+str(modelname)+"_Amin_"+str(Amin)+"_Amax_"+str(Amax)+"_fmin_"+str(fmin)+"_fmax_"+str(fmax)+"kHz_I0_"+str(I0)+"_"+str(Alen)+"_"+str(flen)+".csv"
    Vrest = fhh.get_average_fit(1.0, 0, t_sim=500, t_offset=500, I_0=I0, model=modelname)

    Amplitudes = np.linspace(Amin, Amax, Alen)
    fstim = np.linspace(fmin, fmax, flen)
    R = np.zeros((flen+1,Alen+1), dtype=np.float32)
    R[1:,0] = fstim
    R[0,1:] = Amplitudes

    for idx, f_stim in enumerate(fstim):
        
        def local_orbit(Amp):
            if Amp/(2*np.pi*f_stim) < Alim:
                #
                return fhh.get_average_fit(f_stim, Amp, t_sim=500, t_offset=500, I_0=I0, model=modelname)-Vrest
            else:
                return 0.0

        if __name__ == "__main__":
            print("... Computing orbits")

            pool = mp.Pool(Ntask)
            start = time.time()
            comp_results = pool.map(local_orbit, Amplitudes)
            pool.close()
            end = time.time()
            print("... all simulations for "+str(modelname)+" with fstim = "+str(f_stim)+" finished in %.2f seconds" % (end - start))
            # Buidling matrix to store
            for k, result in enumerate(comp_results):
                R[idx+1, k+1] = result

    np.savetxt(fname, R, delimiter=",")
