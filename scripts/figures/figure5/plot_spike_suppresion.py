import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14
plt.rcParams["mathtext.fontset"] = "stix"

A_min = 10
A_max = 1000
v_min = -100
v_max = 40
A_critical = 497

fname_1 = "bd_Amin_10_Amax_1000_I0_20_fStim_5.0.csv"

fig = plt.figure(figsize=(5, 3))


# load data
Rstrobo = np.loadtxt(fname_1, dtype=np.float32)
ax0 = fig.add_subplot()

span2 = ax0.axvspan(0, A_critical, facecolor='red', alpha=0.1)
t2 = ax0.text(100, 30, 'multiple spikes')

#span3 = ax0.axvspan(A_critical, A_max, facecolor='lime', alpha=0.1)
t3 = ax0.text(620, 30, 'spike suppression')
ax0.plot([A_critical, A_critical], [v_min, v_max], color='red', alpha=0.8, linestyle='--')

ax0.set_xlim([A_min, A_max])
ax0.set_ylim([v_min, v_max])
ax0.grid('both', alpha = 0.2)

# stroboscopic plot itself
# for k in range(Rstrobo.shape[0]):
#     A_stim = Rstrobo[k,0]
#     V_samples = Rstrobo[k,1:]
#     ax0.plot(A_stim*np.ones(len(V_samples)),V_samples,',b')
ax0.plot(Rstrobo[:,0], Rstrobo[:, 1:], ' ,k', c = 'darkblue')
    
ax0.set_ylabel('Membrane potential $v$ ($\mathrm{mV}$)')
ax0.set_xlabel('Stimulation amplitude $A_{\mathrm{stim}}$ ($\mathrm{\mu A/cm^{2}}$)')
plt.subplots_adjust(left=0.15, right=0.96, top=0.99, bottom=0.19)

#plt.tight_layout()
plt.savefig('figure5.png', dpi=400)
plt.show()