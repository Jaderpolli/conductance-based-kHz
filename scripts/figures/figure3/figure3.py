import numpy as np 
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib import cm
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 18
plt.rcParams["mathtext.fontset"] = "stix"


traj1_fname = "traj_IDC_10_fstim_1.0kHz_Amp_75_Length_200_s.csv"
traj2_fname = "traj_IDC_10_fstim_1kHz_Amp_190_Length_200_s.csv"
traj3_fname = "traj_IDC_10_fstim_1kHz_Amp_410_Length_200_s.csv"
traj4_fname = "traj_IDC_10_fstim_1kHz_Amp_600_Length_200_s.csv"

isi_cloud_fname = "isi_m3h_Amin_10_Amax_1000_fStim_1000_I0_10.csv"

v_min = -135
v_max = 60

#######################
## Figure definition ##
#######################
fig = plt.figure(figsize=(14, 8))  # fig, ax = plt.subplots(layout='constrained', figsize=(10,7))
# create grid for different subplots
spec = gridspec.GridSpec(ncols=2, nrows=16,
                         left = 0.07,
                         right = 0.93,
                         top = 0.98,
                         bottom = 0.04,
                         wspace=0.2,
                         hspace=-0.1
)

# all trajectories

ax_traj1 = fig.add_subplot(spec[0:3, 0])
ax_traj2 = fig.add_subplot(spec[4:7, 0])
ax_traj3 = fig.add_subplot(spec[8:11, 0])
ax_traj4 = fig.add_subplot(spec[12:15, 0])

traj1 = np.loadtxt(traj1_fname, delimiter=',', dtype=np.float32)
traj2 = np.loadtxt(traj2_fname, delimiter=',', dtype=np.float32)
traj3 = np.loadtxt(traj3_fname, delimiter=',', dtype=np.float32)
traj4 = np.loadtxt(traj4_fname, delimiter=',', dtype=np.float32)
t = np.linspace(0.001, 200, len(traj1[:, 1]))

imin, imax = 142000, 200000-1
ax_traj1.plot(t[imin:imax], traj1[imin:imax, 0], color='darkblue', lw = 1.2)
ax_traj2.plot(t[imin:imax], traj2[imin:imax, 0], color='darkblue', lw = 1.2)
ax_traj3.plot(t[imin:imax], traj3[imin:imax, 0], color='darkblue', lw = 1.2)
ax_traj4.plot(t[imin:imax], traj4[imin:imax, 0], color='darkblue', lw = 1.2)

ax_traj1.grid('both')
ax_traj2.grid('both')
ax_traj3.grid('both')
ax_traj4.grid('both')
ax_traj1.set_ylim([-85, 60])
ax_traj2.set_ylim([-105, 60])
ax_traj3.set_ylim([-120, 60])
ax_traj4.set_ylim([-135, 60])
ax_traj1.set_xlim([t[imin], t[imax]])
ax_traj2.set_xlim([t[imin], t[imax]])
ax_traj3.set_xlim([t[imin], t[imax]])
ax_traj4.set_xlim([t[imin], t[imax]])
ax_traj3.set_ylabel(r'membrane voltage')
ax_traj2.set_ylabel(r'$v$ (mV)')

ax_traj1.tick_params(labelbottom = False)
ax_traj2.tick_params(labelbottom = False)
ax_traj3.tick_params(labelbottom = False)

ax_traj1.text(160, 34, r'$A_{\rm{stim}} = 75 \, \rm{\mu A/cm^2}$, $A =  11.9 \, \rm{mV}$')
ax_traj2.text(160, 34, r'$A_{\rm{stim}} = 190\, \rm{\mu A/cm^2}$, $A = 30.2 \, \rm{mV}$')
ax_traj3.text(160, 27, r'$A_{\rm{stim}} = 410\, \rm{\mu A/cm^2}$, $A = 65.3 \, \rm{mV}$')
ax_traj4.text(160, 27, r'$A_{\rm{stim}} = 600\, \rm{\mu A/cm^2}$, $A = 95.5 \, \rm{mV}$')
ax_traj4.set_xlabel(r'Time (ms)')


# spike detection
ax_m3h = fig.add_subplot(spec[0:5, 1])
ax_spikes = ax_m3h.twinx()

ax_spikes.plot([160,190], [0,0], ls = 'dashed', color = 'black', alpha = 0.5, label = r'$v = 0 \, \rm{mV}$')
ax_spikes.plot(t[imin:imax], traj4[imin:imax, 0], color='darkblue', lw = 1.0)
ax_spikes.set_ylim([v_min, v_max])
ax_spikes.set_xlim(160, 190)
ax_spikes.legend(loc = 'upper right', frameon = False)
ax_spikes.set_ylabel(r'$v$ (mV)', color='darkblue')
ax_spikes.tick_params(axis='y', labelcolor='darkblue')
ax_m3h.plot([160, 190], [0.12, 0.12], linestyle='dashed', color='red', alpha = 0.5, label = r'$m^3 h =0.12$')
ax_m3h.plot(t[imin:imax], (traj4[imin:imax, 1]**3)*traj4[imin:imax, 2], color='black', linestyle='dashdot')
ax_m3h.set_xlim(160, 190)
ax_m3h.set_ylim([0, 0.22])
#ax_m3h.grid('both')
ax_m3h.legend(loc = 'best', frameon = False)

ax_m3h.set_xlabel(r'Time (ms)')
ax_m3h.set_ylabel(r'$m^3 h = g_{Na}/ \bar{g}_{Na}$')
#ax_m3h.text(170, 0.13, r'$g_{Na}/ \bar{g_{Na}}=0.12$', color='red')

# ISI cloud
ax_ISI = fig.add_subplot(spec[7:15, 1])

R_ISI = np.loadtxt(isi_cloud_fname)
# stroboscopic plot itself
ax_ISI.plot(R_ISI[:,0], R_ISI[:,1:], ' ,k', c = 'darkblue')
ax_ISI.plot([75,75], [0,20], ls = 'dashed', c = 'black', alpha = 0.5)
ax_ISI.plot([190,190], [0,20], ls = 'dashed', c = 'black', alpha = 0.5)
ax_ISI.plot([410,410], [0,20], ls = 'dashed', c = 'black', alpha = 0.5)
ax_ISI.plot([600,600], [0,20], ls = 'dashed', c = 'black', alpha = 0.5)
ax_ISI.set_ylabel(r'ISI (ms)')
ax_ISI.set_xlabel(r'$A_{stim}$ ($\mathrm{\mu A/ cm^2}$)')
ax_ISI.grid('both')
ax_ISI.set_ylim([0.1, 20])
ax_ISI.set_xlim([0,1000])



inset_ISI = inset_axes(ax_ISI,
                    width="30%", # width = 30% of parent_bbox
                    height=2.8, # height : 1 inch
                    loc='upper right')
A_min, A_max = 560, 640
ISI_min, ISI_max = 2, 8

ax_ISI.add_patch(Rectangle((A_min, ISI_min), (A_max-A_min), (ISI_max-ISI_min),color="forestgreen", linewidth=2, linestyle='--', fill = False))


inset_ISI.plot(R_ISI[:,0], R_ISI[:,1:], ',', c = 'darkblue')
inset_ISI.grid('both')
inset_ISI.tick_params(axis='both', which='major', labelsize=10)
inset_ISI.set_xlim([A_min, A_max])
inset_ISI.set_ylim([ISI_min, ISI_max])

xy1_ISI = (A_min, ISI_max)
xy2_ISI = (A_max, ISI_min)
con_zoom1 = ConnectionPatch(xyA=xy1_ISI, xyB=xy1_ISI, coordsA="data", coordsB="data",
                      axesA=ax_ISI, axesB=inset_ISI, color='forestgreen', alpha=0.4, linestyle='--')
con_zoom2 = ConnectionPatch(xyA=xy2_ISI, xyB=xy2_ISI, coordsA="data", coordsB="data",
                      axesA=ax_ISI, axesB=inset_ISI, color='forestgreen', alpha=0.4, linestyle='--')
ax_ISI.add_artist(con_zoom1)
ax_ISI.add_artist(con_zoom2)


# anotations

fig.text(0.1, 0.95, "(a)", fontsize=18)
fig.text(0.1, 0.72, "(b)", fontsize=18)
fig.text(0.1, 0.49, "(c)", fontsize=18)
fig.text(0.1, 0.25, "(d)", fontsize=18)
fig.text(0.72, 0.95, "(e)", fontsize=18)
fig.text(0.55, 0.54, "(f)", fontsize=18)

spec.tight_layout(fig)
plt.savefig("figure3.png", dpi = 200)
plt.show()
