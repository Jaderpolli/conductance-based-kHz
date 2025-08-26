import fhh
import numpy as np
import multiprocessing as mp
import time

# constant A that should lead in the chaotic big branch
A_const = 1500. / (2*np.pi * 1.25)


# from functools import partial

fmin = 0.1
fmax = 1.5
Npamp = 2000
I_0 = 0
reduced = False

length = 1000


Ntask = 7

def local_orbit(f):
    Amp = A_const * 2 * np.pi * f
    print(Amp, '\t\t', f)
    result = fhh.get_orbit(f, Amp, I_0=I_0, reduced=reduced, t_sim=length*(1/f), t_offset=1000*(1/f),)
    #print(Amp, '\t\t', f, '\t\t', result[0:length].shape)
    return result


if __name__ == "__main__":
    Frequencies = np.linspace(fmin, fmax, Npamp)
    print("... Computing orbits")
    fname = "orbits_ISO_A_" + str(A_const) + '.csv'


    pool = mp.Pool(Ntask)
    start = time.time()
    comp_results = pool.map(local_orbit, Frequencies)
    end = time.time()
    print("... all simulations finished in %.2f seconds" % (end - start))
    # Buidling matrix to store
    R = np.zeros((Npamp, length + 1), dtype=np.float32)
    R[:, 0] = Frequencies.transpose()
    for k, result in enumerate(comp_results):
        R[k, 1:] = result[0:length]
    np.savetxt(fname, R, delimiter=",")
