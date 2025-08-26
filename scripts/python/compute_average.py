import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import fhh

A_min = 0
A_max = 3000
Npamp = 300

model = "HH"
I_0 = 0
f_stim = 1.25

Ntask = 4

fname = "average_HH_"+str(A_min)+"_"+str(A_max)+"_fstim_"+str(f_stim)+"_kHz_I0_0.csv"


def local_orbit(Amp):
    print(Amp, I_0)
    return fhh.get_average_fit(f_stim, Amp, t_sim=500, t_offset=500, I_0=I_0, model=model)

if __name__ == "__main__":
    print("... Computing orbits")
    Amplitudes = np.linspace(A_min, A_max, Npamp)


    """
    R = np.zeros(Npamp, dtype=np.float32)
    for k, Amp in enumerate(Amplitudes):
        print(Amp)
        R[k] = 
    """
    pool = mp.Pool(Ntask)
    start = time.time()
    comp_results = pool.map(local_orbit, Amplitudes)
    end = time.time()
    print("... all simulations finished in %.2f seconds" % (end - start))
    # Buidling matrix to store
    R = np.zeros(Npamp, dtype=np.float32)
    for k, result in enumerate(comp_results):
        R[k] = result

    np.savetxt(fname, R, delimiter=",")
