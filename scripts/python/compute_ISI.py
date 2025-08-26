import fhh
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import os

I_0_min_1, I_0_max_1 = -1, 200
dI_1 = 0.1

def get_ISI(I_0):
    t_stop = 1000
    dt = 0.01
    f_stim = 1
    Amp = 0
    t_offset=500
    t, X = fhh.run_forced_RGC2(t_stop, dt, I_0, f_stim, Amp) # change this for other fhh functions

    t = t[(int(t_offset/dt)):]
    t -= t[0]
    X = X[int(t_offset/dt):,:]
    
    v = X[:, 0]
    peaks, _ = sig.find_peaks(v, threshold=-25)
    if len(peaks) == 0:
        mean_ISI = np.NaN
        mean_FS = np.NaN
    else:
        mean_ISI = np.mean(np.diff(t[peaks]))
        mean_FS = 1e3/mean_ISI
    return mean_ISI, mean_FS
    
if __name__ == '__main__':
    all_I0 = np.linspace(I_0_min_1, I_0_max_1, int((I_0_max_1 - I_0_min_1)/dI_1) + 1)  # range(I_0_min, I_0_max, 1)

    results = np.zeros((len(all_I0), 3), dtype=np.float64)

    for k, I_0 in enumerate(all_I0):
        print(I_0)
        results[k, 0] = I_0
        results[k, 1], results[k, 2] = get_ISI(I_0)

    os.makedirs('../../results/raw/mdl_properties', exist_ok=True)
    np.savetxt('../../results/raw/mdl_properties/ISI_RGC2.txt', results) # change filename for each model

    ISI = results[:, 1]
    FS = results[:, 2]


    
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel(r'$I_{stim} = I_{I_0}$')
    ax1.set_ylabel(r'ISI (ms)', color=color)
    ax1.plot(all_I0, ISI, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'$f_{spike}$ (Hz)', color=color)  # we already handled the x-label with ax1
    ax2.plot(all_I0, FS, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    

