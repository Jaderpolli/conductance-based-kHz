import numpy as np
import matplotlib.pyplot as plt
import fhh
from matplotlib import gridspec
from matplotlib import cm
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.rcParams["mathtext.fontset"] = "stix"


v_min = -90
v_max = 30
f_min = 1e1
f_max = 1e4
N_pts = 1e3
N_pts_dec = 10

A_min = 0
A_max = 2500
Npamp = 300

fname_bd_4D = "bd_Amin_10_Amax_1000_I0_0_fStim_1.0.csv"
fname_bd_3D = "h_bd_Amin_0_Amax_700_I0_0_fStim_1.0.csv"
fname_bd_2D = "m_bd_Amin_0_Amax_700_I0_0_fStim_1.0.csv"

fname_avg0 = "average_HH_0_1000_fstim_0.25_kHz_I0_0.csv"
fname_avg1 = "average_HH_0_2500_fstim_1.0_kHz_I0_0.csv"
fname_avg2 = "average_HH_0_1850_fstim_0.5_kHz_I0_0.csv"
fname_avg3 = "average_HH_0_3000_fstim_1.25_kHz_I0_0.csv"

fig = plt.figure(figsize=(14, 8))  # fig, ax = plt.subplots(layout='constrained', figsize=(10,7))
# create grid for different subplots
spec = gridspec.GridSpec(ncols=2, nrows=1,
                         wspace=0.05,
                         hspace=0.05,
                         width_ratios=[3, 2],
)

spec1 = spec[0,0].subgridspec(ncols = 1, nrows = 3, wspace = 0.01)


# compare stroboscopic plots
ax0 = fig.add_subplot(spec1[0, 0])
Rstrobo = np.loadtxt(fname_bd_4D)
# stroboscopic plot itself
# for k in range(Rstrobo.shape[0]):
#     A_stim = Rstrobo[k,0]
#     V_samples = Rstrobo[k,1:]
#     ax0.plot(A_stim*np.ones(len(V_samples)),V_samples,',b')
ax0.plot(Rstrobo[:,0], Rstrobo[:,1:], ' k,', c = 'darkblue')
ax0.set_xlim([0, 700])
ax0.set_ylim([-150, 60])
ax0.set_ylabel('$v$ (mV)')
ax0.set_xticklabels([])
ax0.grid('both')
#ax0.set_title(r'4 dimensions HH')

ax1 = fig.add_subplot(spec1[1, 0])
Rstrobo = np.loadtxt(fname_bd_3D)
# stroboscopic plot itself
# for k in range(Rstrobo.shape[0]):
#     A_stim = Rstrobo[k,0]
#     V_samples = Rstrobo[k,1:]
#     ax1.plot(A_stim*np.ones(len(V_samples)),V_samples,',b')
ax1.plot(Rstrobo[:,0], Rstrobo[:,1:], ' k,', c = 'darkblue')
ax1.set_xlim([0, 700])
ax1.set_ylim([-150, 60])
ax1.set_ylabel('$v$ (mV)')
ax1.set_xticklabels([])
ax1.grid('both')
#ax1.set_title(r'3 dimensions HH')

ax2 = fig.add_subplot(spec1[2, 0])
Rstrobo = np.loadtxt(fname_bd_2D)
# stroboscopic plot itself
# for k in range(Rstrobo.shape[0]):
#     A_stim = Rstrobo[k,0]
#     V_samples = Rstrobo[k,1:]
#     ax2.plot(A_stim*np.ones(len(V_samples)),V_samples,',b')
ax2.plot(Rstrobo[:,0], Rstrobo[:,1:], ' k,', c = 'darkblue')
ax2.set_xlim([0, 700])
ax2.set_ylim([-150, 60])
ax2.set_ylabel('$v$ (mV)')
ax2.set_xlabel('Stimulation amplitude $A_{stim}$ ($\mathrm{\mu A/cm^2}$)')
ax2.grid('both')
#ax2.set_title(r'2 dimensions HH')

spec2 = spec[0,1].subgridspec(ncols = 1, nrows = 2, hspace = 0.3)

# plot average values
ax3 = fig.add_subplot(spec2[1, 0])
R_avg0 = np.loadtxt(fname_avg0)
Amplitudes0 = np.linspace(A_min, 1000, 200)
R_avg1 = np.loadtxt(fname_avg1)
Amplitudes1 = np.linspace(A_min, A_max, Npamp)
R_avg2 = np.loadtxt(fname_avg2)
Amplitudes2 = np.linspace(A_min, 1850, 200)
R_avg3 = np.loadtxt(fname_avg3)
Amplitudes3 = np.linspace(A_min, 3000, 300)
ax3.plot(R_avg0, Amplitudes0, color="gold", label=r"0.25 kHz")
ax3.plot(R_avg2, Amplitudes2, color="darkorange", label=r"0.5 kHz")
ax3.plot(R_avg1, Amplitudes1, color="firebrick", label=r"1 kHz")
ax3.plot(R_avg3, Amplitudes3, color="darkorchid", label=r"1.25 kHz")
ax3.grid('both')
ax3.set_xlabel(r" $\langle v \rangle$ (mV)")
ax3.set_ylabel(r"$A_{stim}$ ($\mathrm{\mu A/cm^2}$)")
ax3.set_ylim([0, 3000])
ax3.set_xlim([v_min, v_max])
ax3.legend(loc='lower right', frameon = False)
ax3.plot([-65,-65],[0,3000], linestyle = "dashed", color = 'gray', alpha = 0.5)
#span1 = ax3.axvspan(-90,-65, facecolor='blue', alpha=0.1)
#span2 = ax3.axvspan(-65,30, facecolor = 'red', alpha = 0.1)

# plot particules frequency versus voltage evolution
v_scale = np.linspace(v_min, v_max, int(N_pts))
f_scale = np.logspace(np.log10(f_min), np.log10(f_max), int(np.log10(f_max/f_min)*N_pts_dec)) 
w_scale = 2 * np.pi * f_scale

tau_m = fhh.tau_m(v_scale)
tau_h = fhh.tau_h(v_scale)
tau_n = fhh.tau_n(v_scale)

f_m = 1./(2 * np.pi * tau_m)
f_h = 1./(2 * np.pi * tau_h)
f_n = 1./(2 * np.pi * tau_n)

ax4 = fig.add_subplot(spec2[0, 0])

#ax4.plot([v_min, v_max], [0.25, 0.25], linestyle="dashed", color="gold", alpha=0.5)
#ax4.text(-100, 0.25, r'$f_{stim} = 0.25kHz$', fontsize=12, color='gray')
#ax4.plot([v_min, v_max], [0.5, 0.5], linestyle="dashed", color="darkorange", alpha=0.5)
#ax4.text(-100, 0.5, r'$f_{stim} = 0.5kHz$', fontsize=12, color='gray')
#ax4.plot([v_min, v_max], [1.0, 1.0], linestyle="dashed", color="firebrick", alpha=0.5)
#ax4.text(-100, 1.0, r'$f_{stim} = 1.0kHz$', fontsize=12, color='gray')
#ax4.plot([v_min, v_max], [1.25, 1.25], linestyle="dashed", color="darkorchid", alpha=0.5)
#ax4.text(-70, 1.25, r'$f_{stim} = 1.25kHz$', fontsize=12, color='gray')


ax4.plot(v_scale, f_m, label=r'$f_m$', color='forestgreen')
ax4.plot(v_scale, f_h, label=r'$f_h$', color='palegreen')
ax4.plot(v_scale, f_n, label=r'$f_n$', color='cornflowerblue')

fm_min = np.min(f_m)
print('fm minimum value ', fm_min)
fh_min = np.min(f_h)
print('fh minimum value ', fh_min)
v_fm_min = v_scale[np.argmin(f_m)]
v_fh_min = v_scale[np.argmin(f_h)]

#ax4.scatter(v_fm_min, fm_min, color='forestgreen', alpha = 0.5)
#ax4.plot([v_min, v_max], [fm_min, fm_min], color='forestgreen', linestyle='dotted', alpha=0.5)
#ax3.plot([v_fm_min, v_fm_min], [0, 2000], color='forestgreen', linestyle='dotted', alpha=0.5)
#ax3.text(v_fm_min, 1500, r'$v|_{f_{m} = f_{m\, min}}$', rotation=90, color='gray')
#ax4.plot([v_fm_min, v_fm_min], [0.01, 5], color='forestgreen', linestyle='dotted', alpha=0.5)


#ax4.scatter(v_fh_min, fh_min, color='palegreen', alpha = 0.5)
#ax4.plot([v_min, v_max], [fh_min, fh_min], color='palegreen', linestyle='dotted', alpha=0.5)
#ax3.plot([v_fh_min, v_fh_min], [0, 2000], color='palegreen', linestyle='dotted', alpha=0.5)
#ax3.text(v_fh_min, 300, r'$v|_{f_{h} = f_{h\, min}}$', rotation=90, color='gray')
#ax4.plot([v_fh_min, v_fh_min], [0.01, 5], color='palegreen', linestyle='dotted', alpha=0.5)


ax4.grid('both')
ax4.set_yscale('log')
ax4.set_xlim([v_min, v_max])
ax4.set_ylim([0.005, 5])
ax4.set_xlabel(r'$v$ (mV)')
ax4.set_ylabel(r'$f_x$ (kHz)')
ax4.legend(loc='lower right', frameon = False)



#xy0_4 = (v_fm_min, 1e-2)
#xy0_3 = (v_fm_min, 2000)
#con1 = ConnectionPatch(xyA=xy0_4, xyB=xy0_3, coordsA="data", coordsB="data",
#                      axesA=ax4, axesB=ax3, color="forestgreen", linestyle='dotted', alpha = 0.4)
#ax4.add_artist(con1)

#xy1_4 = (v_fh_min, 1e-2)
#xy1_3 = (v_fh_min, 2000)
#con2 = ConnectionPatch(xyA=xy1_4, xyB=xy1_3, coordsA="data", coordsB="data",
#                      axesA=ax4, axesB=ax3, color="palegreen", linestyle='dotted', alpha = 0.4)
#ax4.add_artist(con2)


#### PAPER LABELS ##
fig.text(0.01, 0.95, "(a)", fontsize=20)
fig.text(0.01, 0.63, "(b)", fontsize=20)
fig.text(0.01, 0.31, "(c)", fontsize=20)
fig.text((3./5.)-0.01, 0.95, "(d)", fontsize=20)
fig.text((3./5.)-0.01, 0.48, "(e)", fontsize=20)


spec.tight_layout(fig)
plt.savefig('figure7.png', dpi=500)
plt.show()