import numpy as np
import matplotlib.pyplot as plt


ISI_file_1 = "../../../../results/raw/mdl_properties/ISI_RGC1.txt"
ISI_file_2 = "../../../../results/raw/mdl_properties/ISI_RGC2.txt"

data = np.loadtxt(ISI_file_1)
all_I0_RGC1 = data[:, 0]
ISI_RGC1 = data[:, 1]
FS_RGC1 = data[:, 2]

I_min_RGC1 = all_I0_RGC1[0]
I_max_RGC1 = all_I0_RGC1[-1]
ISI_min_RGC1 = np.nanmin(ISI_RGC1)
ISI_max_RGC1 = np.nanmax(ISI_RGC1)
f_max_RGC1 = np.nanmax(FS_RGC1)*1.2

data = np.loadtxt(ISI_file_2)
all_I0_RGC2 = data[:, 0]
ISI_RGC2 = data[:, 1]
FS_RGC2 = data[:, 2]

I_min_RGC2 = all_I0_RGC2[0]
I_max_RGC2 = all_I0_RGC2[-1]
ISI_min_RGC2 = np.nanmin(ISI_RGC2)
ISI_max_RGC2 = np.nanmax(ISI_RGC2)
f_max_RGC2 = np.nanmax(FS_RGC2)*1.2

fig, ax = plt.subplots(1, 2, figsize=(9, 4))

color = 'grey'
ax[0].set_ylabel(r'ISI (ms)', color=color)
ax[0].set_xlim([-5, 1.2*I_max_RGC1])
ax[0].set_ylim([ISI_min_RGC1, ISI_max_RGC1])
ax[0].plot(all_I0_RGC1, ISI_RGC1, color=color)
ax[0].tick_params(axis='y', labelcolor=color)
ax[0].grid()
ax[0].set_title("RGC 1 model")
ax[0].set_xlabel(r'$I_{stim} = I_{0}$ ($\mu$A/cm$\,^2$)')

ax01 = ax[0].twinx()  # instantiate a second Axes that shares the same x-axis
color = 'black'
ax01.set_ylabel(r'$f_{osc}$ (Hz)', color=color)  # we already handled the x-label with ax1
ax01.plot(all_I0_RGC1, FS_RGC1, color=color, linestyle='dashed')
ax01.set_ylim([0, f_max_RGC1])
ax01.tick_params(axis='y', labelcolor=color)



color = 'grey'
ax[1].set_ylabel(r'ISI (ms)', color=color)
ax[1].set_xlim([-20, 1.2*I_max_RGC2])
ax[1].set_ylim([ISI_min_RGC2, ISI_max_RGC2])
ax[1].plot(all_I0_RGC2, ISI_RGC2, color=color)
ax[1].tick_params(axis='y', labelcolor=color)
ax[1].grid()
ax[1].set_title("RGC 2 model")
ax[1].set_xlabel(r'$I_{stim} = I_{0}$ ($\mu$A/cm$\,^2$)')

ax11 = ax[1].twinx()  # instantiate a second Axes that shares the same x-axis
color = 'black'
ax11.set_ylabel(r'$f_{osc}$ (Hz)', color=color)  # we already handled the x-label with ax1
ax11.plot(all_I0_RGC2, FS_RGC2, color=color, linestyle='dashed')
ax11.set_ylim([0, f_max_RGC2])
ax11.tick_params(axis='y', labelcolor=color)

plt.tight_layout()

fig.savefig('16_RGC.pdf')
plt.show()