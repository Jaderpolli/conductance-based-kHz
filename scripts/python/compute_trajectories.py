import fhh
import numpy as np
import matplotlib.pyplot as plt

transient = 0
t_stop = 200
dt = 0.001
N = int(t_stop / dt)
IDC = [10]#, 5, 10]
Fstims = [1.0]
Amplitudes = [75]

trec = 50

Nrec = int(trec / dt)

for I_0 in IDC:
    for f_stim in Fstims:
        for Amp in Amplitudes:
            print('Amplitude : ' + str(Amp) + 'I_0 = ' + str(I_0))
            t, X = fhh.run_forced_HH(
                transient + t_stop,
                dt,
                I_0,
                f_stim,
                Amp,
            )
            # get rid of transient
            t = t[(int(transient/dt)):]
            t -= t[0]
            X = X[int(transient/dt):,:]
            fig, ax = plt.subplots()
            ax.plot(t,X[:,0], label='forced')
            plt.show()
            fname = 'traj_IDC_'+str(I_0)+'_fstim_'+str(f_stim)+'kHz_Amp_'+str(Amp)+'_Length_'+str(t_stop)+'_s.csv'
            np.savetxt(fname, X, delimiter=",")


