import fhh
import numpy as np
import multiprocessing as mp
import time

# from functools import partial

f_stim = .5
Amins = [0]  # , 570, 620]
Amaxs = [1800]  # , 650, 621]
Npamp = 3601
I_0 = 0
reduced = False
which = "3D" # either '2D', '3D', 'fast'

Ntask = 12

def local_orbit(Amp):
    print(Amp, I_0)
    return fhh.get_orbit(f_stim, Amp, I_0=I_0, reduced=reduced, which=which)


if __name__ == "__main__":
    print("... Computing orbits")
    for k in range(len(Amins)):
        Amin = Amins[k]
        Amax = Amaxs[k]
        Amplitudes = np.linspace(Amin, Amax, Npamp)

        if reduced == True:
            if which == "fast":
                fname = (
                    "orbits_Amin_"
                    + str(Amin)
                    + "_Amax_"
                    + str(Amax)
                    + "_I0_"
                    + str(I_0)
                    + "_reducedHH_fast.csv"
                )
            elif which == "3D":
                fname = (
                    "orbits_Amin_"
                    + str(Amin)
                    + "_Amax_"
                    + str(Amax)
                    + "_I0_"
                    + str(I_0)
                    + "_reducedHH_3D.csv"
                )
            else:
                fname = (
                    "orbits_Amin_"
                    + str(Amin)
                    + "_Amax_"
                    + str(Amax)
                    + "_I0_"
                    + str(I_0)
                    + "_reducedHH.csv"
                )
        else:
            fname = (
                "orbits_Amin_"
                + str(Amin)
                + "_Amax_"
                + str(Amax)
                + "_I0_"
                + str(I_0)
                + "_fstim_"
                #+ '500Hz'
                + str(f_stim)
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
