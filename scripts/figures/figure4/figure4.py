import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.gridspec import GridSpec

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams["mathtext.fontset"] = "stix"

# Arquivos
fname_1 = "bd_Amin_10_Amax_3000_I0_0_fStim_1.25.csv"
fname_2 = "lyapunov_Amin_10_Amax_3000_fStim_1250_I0_0_size_5000.csv"
fname_3 = "bd_Amin_1850_Amax_2200_I0_0_fStim_1.25.csv"
fname_4 = "spectra_I0_0_fstim_1250Hz_0_1300.csv"
fscale_name = "fscale_1250Hz.csv"

# Parâmetros
A_min, A_max = 0, 3000
v_min, v_max = -160, 40
v_zone_min, v_zone_max = -150, 25
L_min, L_max = -1.5, 0.5
A_min_zoom, A_max_zoom = 1930, 2125
v_min_zoom, v_max_zoom = -110, -50

# Dados
Rstrobo = np.loadtxt(fname_1, dtype=np.float32)
L_values = np.loadtxt(fname_2, dtype=np.float32)
R = np.loadtxt(fname_4, dtype=np.float32, delimiter=",")
f = np.loadtxt(fscale_name, dtype=np.float32, delimiter=",")
amplitudes = R[:, 0]
surface = R[:, 1:]

# Figura com layout mais limpo
fig = plt.figure(figsize=(11, 8), constrained_layout=True)
gs = GridSpec(5, 4, width_ratios=[20, 1.2, 1, 10], height_ratios=[2,0.1,1, 1, 1])

#####################
 (a)
#####################
ax0 = fig.add_subplot(gs[0:3, 0])
ax0.set_xlim([A_min, A_max])
ax0.set_ylim([v_min, v_max+10])
ax0.grid(True)
ax0.axvspan(0, 140.5, facecolor='lavender', alpha=0.35)
ax0.axvspan(140.5, A_max, facecolor='red', alpha=0.12)
ax0.plot(Rstrobo[:,0], Rstrobo[:,1:], ',k', c='darkblue')
ax0.add_patch(Rectangle((A_min_zoom, v_min_zoom), A_max_zoom - A_min_zoom, v_max_zoom - v_min_zoom,
                        edgecolor='forestgreen', linestyle='--', fill=False, linewidth=1.5))
# Pn labels (reduzido)
ax0.plot([2815., 2815], [v_zone_min, v_zone_max], color='gray', linestyle='dotted')
ax0.text(2840, 19, r'$P_2$')

ax0.plot([2120, 2120], [v_zone_min, v_zone_max], color='gray', linestyle='dotted')
ax0.plot([2348, 2348], [v_zone_min, v_zone_max], color='gray', linestyle='dotted')
ax0.text(2160, 19, r'$P_3$')

ax0.plot([830, 830], [v_zone_min, v_zone_max], color='gray', linestyle='dotted')
ax0.plot([1168, 1168], [v_zone_min, v_zone_max], color='gray', linestyle='dotted')
ax0.text(950, 19, r'$P_4$')


ax0.text(340, 19, r'$\dots P_5$')
ax0.set_ylabel('$v$ (mV)')
ax0.tick_params(labelbottom=False)
ax0.annotate("(a)", xy=(0.02, 0.95), xycoords="axes fraction", ha="left")

#####################
# (b)
#####################
ax1 = fig.add_subplot(gs[3, 0])
ax1.set_xlim([A_min, A_max])
ax1.set_ylim([L_min, L_max])
ax1.plot([A_min, A_max], [0, 0], color='k', linewidth=1, alpha=0.5)
ax1.plot(L_values[:, 0], L_values[:, 1], color='darkblue', linewidth = 0.75)
ax1.grid(True)
ax1.set_ylabel('$\lambda$')
#ax1.set_xlabel(r'$A_{\rm{stim}}$ ($\rm{\mu A/cm^2}$)')
ax1.tick_params(labelbottom=False)
ax1.annotate("(b)", xy=(0.02, 0.82), xycoords="axes fraction", ha="left")

#####################
# (d)
#####################
ax2 = fig.add_subplot(gs[4, 0])
pcm = ax2.pcolormesh(amplitudes, f, np.log10(surface).T, cmap=cm.nipy_spectral, rasterized=True)
ax2.set_ylabel(r"$f$ (kHz)")
ax2.set_xlabel(r"$A_{\rm{stim}}$ ($\mu$A/cm$^2$)")
ax2.annotate("(c)",xy=(0.02, 0.82), xycoords="axes fraction",ha="left")

#####################
(c)
#####################
cbar_ax = fig.add_subplot(gs[4, 1])
cbar = fig.colorbar(pcm, cax=cbar_ax)
cbar.ax.text(-0.5, 0.5, r"$\log_{10} [PSD]$", rotation=90,
             va='center', ha='center', transform=cbar.ax.transAxes)
cbar.set_ticks([-10, -5, 0, 5])
#cbar.set_ticklabels([r"$10^{-10}$", r"$10^{-5}$", r"$10^0$", r"$10^5$"])
#cbar.ax.set_ylabel(r"PSD $\left(\rm{mV^2/Hz} \right)$")
#cbar.ax.set_label(r"PSD ($\left(\rm{mV^2/Hz} \right)$)")

#############################
Rstrobo_zoom = np.loadtxt(fname_3, dtype=np.float32, delimiter='\t')
ax3 = fig.add_subplot(gs[0:1, 2:4])
ax3.plot(Rstrobo_zoom[:,0], Rstrobo_zoom[:,1:], ',k', c='darkblue')
ax3.set_xlim([A_min_zoom, A_max_zoom])
ax3.set_xticks([A_min_zoom, A_max_zoom])
ax3.set_ylim([v_min_zoom, v_max_zoom])
ax3.set_yticks([v_min_zoom, v_max_zoom])
ax3.grid(True)
ax3.annotate("(d)", xy=(0.5, 0.9), xycoords="axes fraction", ha="left")

# Conexões com ax0
for y in [v_min_zoom, v_max_zoom]:
    con = ConnectionPatch(
        xyA=(A_max_zoom, y), xyB=(A_min_zoom, y),
        coordsA="data", coordsB="data",
        axesA=ax0, axesB=ax3,
        color="forestgreen", linestyle='--', alpha=0.4
    )
    ax3.add_artist(con)

#############################
# (e), (f), (g)
#############################
axes_spectra = [fig.add_subplot(gs[i, 3]) for i in [2, 3, 4]]
target_As = [1000, 480, 2500]
labs = ["(e)", "(f)", "(g)"]

for i, (A_target, ax) in enumerate(zip(target_As, axes_spectra)):
    idx = np.argmin(np.abs(amplitudes - A_target))
    spectrum = surface[idx, :]
    ax.plot(f, np.log10(spectrum), color='darkblue', linewidth = 0.75)
    ax.set_xlim(f[0], f[-1])
    ax.set_ylim(min(np.log10(spectrum)),max(np.log10(spectrum))+2.5)
    ax.grid(True)
    ax.annotate(f"$A_{{stim}}$ = {amplitudes[idx]:.0f} $\\mu$A/cm$^2$",
                xy=(0.98, 0.825), xycoords='axes fraction',
                ha='right', fontsize=14)
    ax.annotate(labs[i], xy=(0.02, 0.825), xycoords="axes fraction", ha="left")
    if i < 2:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(r"$f$  (kHz)")
    if i == 1:
        ax.set_ylabel(r"$\log_{10} \, [PSD] \, (\rm{mV^2 Hz^{-1}})$")

# Export
plt.savefig("figure4.png", dpi=500)
plt.show()
