import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import cm
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.rcParams["mathtext.fontset"] = "stix"


v_max = 20
v_min = -150

fmin = 0.05
fmax = 1.5

fm_min = 0.3174005219037049

f_spike_max = 0.1690273093841642265

fname_1 = "orbits_ISO_A_190.9859317102744.csv"


# load data
Rstrobo = np.loadtxt(fname_1, delimiter=',', dtype=np.float32)


plt.figure(figsize=(7,5))
#plt.axvspan(f_spike_max, fm_min, facecolor='gray', alpha=0.1)
#plt.axvspan(fm_min, fmax, facecolor='forestgreen', alpha=0.1)
#plt.text(0.4, 10, r'$f_{stim} > f_{m\, min}$')
#plt.text(0.2, -80, r'$f_{spike\, max} <f_{stim} < f_{m\, min}$', rotation=90)
# for k in range(Rstrobo.shape[0]):
#     F_stim = Rstrobo[k,0]
#     V_samples = Rstrobo[k,1:]
#     plt.plot(F_stim*np.ones(len(V_samples)),V_samples,',b')
plt.plot(Rstrobo[:,0], Rstrobo[:,1:], ' k,', c = 'darkblue')
plt.grid('both')
plt.plot([fm_min, fm_min], [v_min, v_max], linestyle='dashed', color='forestgreen', alpha=0.5)
plt.plot([f_spike_max, f_spike_max], [v_min, v_max], linestyle='dashed', color='blue', alpha=0.5)
plt.xlabel(r'$f_{\mathrm{stim}}$ (kHz)')
plt.ylabel(r'$v$ (mV)')
plt.ylim([v_min, v_max])
plt.xlim([0.05, 1.5])

plt.tight_layout()

plt.savefig('figure8.png', dpi=400)
#plt.show()
