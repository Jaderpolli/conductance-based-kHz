import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 12
plt.rc('legend', fontsize = 10)
plt.rcParams["mathtext.fontset"] = "stix"

fig, ax = plt.subplots(2,1, figsize = (4,5))
data = np.loadtxt("ISI_HH4D.txt")

all_I0 = data[:, 0]
ISI = data[:, 1]
FS = data[:, 2]

I_0_min = 0
I_0_max = all_I0[-1]
ISI_min = np.nanmin(ISI)
ISI_max = np.nanmax(ISI)
f_max = np.nanmax(FS)*1.2

color1 = 'grey'

p1 = ax[0].plot(all_I0, ISI, color=color1)

color = 'grey'
ax[0].set_ylabel(r'ISI (ms)', color=color1)
ax[0].set_xlim([I_0_min, I_0_max])
ax[0].set_ylim([ISI_min, ISI_max])

ax[0].tick_params(axis='y', labelcolor=color1)

ax2 = ax[0].twinx()  # instantiate a second Axes that shares the same x-axis

color2 = 'black'
ax2.set_ylabel(r'$f_{\rm{osc}}$ (Hz)', color=color2)  # we already handled the x-label with ax1
ax2.plot(all_I0, FS, color=color2, linestyle='dashed')
ax2.set_ylim([0, f_max])
ax2.tick_params(axis='y', labelcolor=color2)

print('fmax')
print(np.nanmax(FS))


stable = np.loadtxt("./stable.dat", delimiter=",")

I_01 = stable[0:29, 1]
V1 = stable[0:29, 2]

I_02 = stable[30:, 1]
V2 = stable[30:, 2]

p2 = ax[1].plot(I_01, V1, label = 'Stable Eq.', c = 'red')

ax[1].plot(I_02, V2, c = 'red')


unstable = np.loadtxt("./unstable.dat", delimiter=",")

I_0 = unstable[:, 1]
V = unstable[:, 2]

ax[1].plot(I_0, V, label = 'Unstable Eq.', c = 'black')

cyclestable = np.loadtxt("./cycle-stable.dat", delimiter=",")

I_0 = cyclestable[:, 1]
V1 = cyclestable[:, 2]
V2 = cyclestable[:, 3]

ax[1].plot(I_0, V1, label = 'Stable LC', mfc = 'none', c = 'y', marker = 'o', ls = "", ms = 2)
ax[1].plot(I_0, V2, mfc = 'none', c = 'y', marker = 'o', ls = "", ms = 2)

cycleunstable = np.loadtxt("./cycle-unstable.dat", delimiter=",")

I_0 = cycleunstable[:, 1]
V1 = cycleunstable[:, 2]
V2 = cycleunstable[:, 3]

ax[1].plot(I_0, V1, label = 'Untable LC',mfc = 'none', c = 'b', marker = 'o', ls = "", ms = 2)
ax[1].plot(I_0, V2, mfc = 'none', c = 'b', marker = 'o', ls = "", ms = 2)

ax[1].set_xlim([0, np.nanmax(stable[:, 1])])
ax[0].set_xlim([0, np.nanmax(stable[:, 1])])
ax[1].set_ylim([np.nanmin(cycleunstable[:, 3])-10, np.nanmax(cyclestable[:, 2])+10])
ax[1].set_ylabel(r'$v$')
ax[1].set_xlabel(r'$I_0$')


ax[1].text(12,-20, 'Subcritical \nHopf Bifurcation', size = 10, multialignment = 'center')
ax[1].annotate(r'$I_{H_1} = 9.775357$', xytext = (16,-30), xy = (10, -59), size = 10, arrowprops = dict(arrowstyle = '->'))

ax[1].text(123,-82, 'Supercritical \nHopf Bifurcation', size = 10, multialignment = 'center')
ax[1].annotate(r'$I_{H_2} = 154.5226$', xytext = (130,-65), xy = (155, -42), size = 10, arrowprops = dict(arrowstyle = '->'))

ax[0].text(80,18, '(c)')
ax[1].text(80,30, '(d)')

ax[1].legend(loc = 'best', frameon = False)
ax[1].grid('both', alpha = 0.3)
ax[0].grid('both', alpha = 0.3)

fig.subplots_adjust(top=0.985,
                    bottom=0.096,
                    left=0.155,
                    right=0.84,
                    hspace=0.162,
                    wspace=0.19)
#fig.tight_layout()

plt.savefig("figure2.pdf")
plt.show()