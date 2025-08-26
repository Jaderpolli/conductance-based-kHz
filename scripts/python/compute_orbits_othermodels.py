import fhh
import numpy as np
import multiprocessing as mp
import time

# from functools import partial

f_stim = .2
Amin = 0
Amax = 100
Npamp = 5001
I_0 = 0
model = "TermanGPe" # either '2D', '3D', 'fast

Ntask = 14

def local_orbit(Amp):
    print(Amp, I_0)
    return fhh.get_orbit_othermodels(f_stim, Amp, I_0=I_0, model=model)


if __name__ == "__main__":
    print("... Computing orbits")
    Amplitudes = np.linspace(Amin, Amax, Npamp)

    fname = (
        "orbits_Amin_"
        + str(Amin)
        + "_Amax_"
        + str(Amax)
        + "_I0_"
        + str(I_0)
        + "_"
        + model
        + ".csv"
    )

    #print(Amin, Amax, Amplitudes)
    pool = mp.Pool(Ntask)
    start = time.time()
    comp_results = pool.map(local_orbit, Amplitudes)
    end = time.time()
    print("... all simulations finished in %.2f seconds" % (end - start))
    # Buidling matrix to store
    R = np.zeros((Npamp, len(comp_results[0]) + 1), dtype=np.float32)
    R[:, 0] = Amplitudes.transpose()
    for k, result in enumerate(comp_results):
        R[k, 1:] = result
    np.savetxt(fname, R, delimiter=",")
