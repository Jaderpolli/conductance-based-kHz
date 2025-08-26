
import fhh
import numpy as np
#import matplotlib.pyplot as plt
import multiprocessing as mp
import time

# processing
t_stop = 10000
dt = 1.e-3
I_0 = 0
f_stim = 1.25
f_max = 1.35
Amax = 3000
Npamp = 3000
Ntask = 5



Amplitudes = np.linspace(1, Amax, Npamp)

path_offset = ''

fspectra_name = path_offset + 'spectra_I0_0_fstim_1250Hz_0_1300.csv'
fscale_name = path_offset + 'fscale_1250Hz.csv'

def get_spectrum(Amp, t_offset=500):
    print(Amp)
    # compute output
    t, X = fhh.run_forced_HH(t_stop+t_offset,
        dt,
        I_0,
        f_stim,
        Amp
    )
    # removing transient
    t = t[(int(t_offset/dt)):]
    t -= t[0]
    X = X[int(t_offset/dt):,:]
    # compute spectrum
    f, Xf = fhh.get_spectral_vector(t, X)

    good_args = np.where(f<f_max)
    return Xf[good_args,0][0]

def get_fscale(t_offset=500):
    # compute output
    t, X = fhh.run_forced_HH( t_stop+t_offset,
        dt,
        I_0,
        f_stim,
        0
    ) # removing transient
    t = t[(int(t_offset/dt)):]
    t -= t[0]
    X = X[int(t_offset/dt):,:]
    # compute spectrum
    f, Xf = fhh.get_spectral_vector(t, X)

    good_args = np.where(f<f_max)
    return f[good_args]

if __name__ == '__main__':
    
    # computing freq
    print('... Computing frequency vector')
    f = get_fscale()
    np.savetxt(fscale_name, f, delimiter=',')
    print('... done')
    
    ##
    print('... Computing spectra')
    pool = mp.Pool(Ntask)
    start = time.time()
    comp_results = pool.map(get_spectrum, Amplitudes)
    end = time.time()
    print('... all simulations finished in %.2f seconds'%(end-start))
    # Building matrix to store
    R = np.zeros((Npamp, len(comp_results[0])+1), dtype=np.float32)
    print(R.shape)
    R[:, 0] = Amplitudes.transpose()
    for k, result in enumerate(comp_results):
        R[k, 1:] = result
    np.savetxt(fspectra_name, R, delimiter=",")
