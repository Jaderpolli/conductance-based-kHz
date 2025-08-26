import numpy as np
import matplotlib.pyplot as plt

bifurcation_file = "./02_hh3Dm.dat"
ISI_file = "../../../../results/raw/mdl_properties/ISI_HH3Dm.txt"


def read_xppdat(filename_bifurcation, filename_ISI):
    STABLE = 1
    UNSTABLE = 2
    STABLE_LC = 3
    UNSTABLE_LC = 4

    CURVES = []

    I_stim = np.loadtxt(
        filename_bifurcation, delimiter=" ", usecols=0, dtype=np.float64
    )
    V_eq1 = np.loadtxt(filename_bifurcation, delimiter=" ", usecols=1, dtype=np.float64)
    V_eq2 = np.loadtxt(filename_bifurcation, delimiter=" ", usecols=2, dtype=np.float64)
    Eq_type1 = np.loadtxt(
        filename_bifurcation, delimiter=" ", usecols=3, dtype=np.int32
    )
    Eq_type2 = np.loadtxt(
        filename_bifurcation, delimiter=" ", usecols=3, dtype=np.int32
    )

    mask_stable_1 = np.argwhere(Eq_type1 == STABLE)
    mask_stable_2 = np.argwhere(Eq_type2 == STABLE)
    mask_unstable_1 = np.argwhere(Eq_type1 == UNSTABLE)
    mask_unstable_2 = np.argwhere(Eq_type2 == UNSTABLE)
    mask_unstable_1 = np.argwhere(Eq_type1 == UNSTABLE)
    mask_unstable_2 = np.argwhere(Eq_type2 == UNSTABLE)
    mask_stable_lc_1 = np.argwhere(Eq_type1 == STABLE_LC)
    mask_stable_lc_2 = np.argwhere(Eq_type2 == STABLE_LC)
    mask_unstable_lc_1 = np.argwhere(Eq_type1 == UNSTABLE_LC)
    mask_unstable_lc_2 = np.argwhere(Eq_type2 == UNSTABLE_LC)

    stable_Istim = np.concatenate((I_stim[mask_stable_1], I_stim[mask_stable_2]))
    stable_V_eq = np.concatenate((V_eq1[mask_stable_1], V_eq2[mask_stable_2]))
    stable_sorter = np.argsort(stable_V_eq, axis=0)
    stable_Istim = stable_Istim[stable_sorter][:, :, 0]
    stable_V_eq = stable_V_eq[stable_sorter][:, :, 0]

    unstable_Istim = np.concatenate((I_stim[mask_unstable_1], I_stim[mask_unstable_2]))
    unstable_V_eq = np.concatenate((V_eq1[mask_unstable_1], V_eq2[mask_unstable_2]))
    unstable_sorter = np.argsort(unstable_V_eq, axis=0)
    unstable_Istim = unstable_Istim[unstable_sorter][:, :, 0]
    unstable_V_eq = unstable_V_eq[unstable_sorter][:, :, 0]

    unstable_Imin = 7.74 # np.min(unstable_Istim)
    unstable_Imax = np.max(unstable_Istim)
    where_trully_unstable = np.where(unstable_Istim >= unstable_Imin)
    unstable_Istim = unstable_Istim[where_trully_unstable]
    unstable_V_eq = unstable_V_eq[where_trully_unstable]
    where_low = np.where(stable_Istim <= unstable_Imin)
    where_more = np.where(stable_Istim >= unstable_Imin)
    curve_stable_1_Istim = stable_Istim[where_low]
    curve_stable_1_V_eq = stable_V_eq[where_low]
    curve_stable_2_Istim = stable_Istim[where_more]
    curve_stable_2_V_eq = stable_V_eq[where_more]

    stable_lc_Istim = np.concatenate(
        (I_stim[mask_stable_lc_1], I_stim[mask_stable_lc_2])
    )
    stable_lc_V_eq = np.concatenate((V_eq1[mask_stable_lc_1], V_eq2[mask_stable_lc_2]))

    unstable_lc_Istim = np.concatenate(
        (I_stim[mask_unstable_lc_1], I_stim[mask_unstable_lc_2])
    )
    unstable_lc_V_eq = np.concatenate(
        (V_eq1[mask_unstable_lc_1], V_eq2[mask_unstable_lc_2])
    )

    data = np.loadtxt(filename_ISI)
    all_I0 = data[:, 0]
    ISI = data[:, 1]
    FS = data[:, 2]

    ISI_min = np.nanmin(ISI)
    ISI_max = np.nanmax(ISI)
    f_max = np.nanmax(FS)*1.2
    print(ISI_min)
    print(ISI_max)


    Imin = -10
    Imax = np.max(I_stim)

    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    color = 'grey'
    ax[0].set_ylabel(r'ISI (ms)', color=color)
    ax[0].set_xlim([Imin, Imax])
    ax[0].set_ylim([ISI_min, ISI_max])
    ax[0].plot(all_I0, ISI, color=color)
    ax[0].tick_params(axis='y', labelcolor=color)
    ax[0].grid()
    fig.suptitle(r"3D-$m$-reduced HH")

    ax2 = ax[0].twinx()  # instantiate a second Axes that shares the same x-axis
    color = 'black'
    ax2.set_ylabel(r'$f_{\mathrm{osc}}$ (Hz)', color=color)  # we already handled the x-label with ax1
    ax2.plot(all_I0, FS, color=color, linestyle='dashed')
    ax2.set_ylim([0, f_max])
    ax2.tick_params(axis='y', labelcolor=color)
    
    
    ax[1].plot(curve_stable_1_Istim, curve_stable_1_V_eq, label="Stable Eq.", c="red")
    ax[1].plot(curve_stable_2_Istim, curve_stable_2_V_eq, c="red")
    ax[1].plot(unstable_Istim, unstable_V_eq, label="Unstable Eq.", c="black")
    ax[1].plot(
        stable_lc_Istim,
        stable_lc_V_eq,
        label="Stable LC",
        mfc="none",
        c="y",
        marker="o",
        ls="",
        ms=2,
    )
    ax[1].plot(
        unstable_lc_Istim,
        unstable_lc_V_eq,
        label="Unstable LC",
        mfc="none",
        c="b",
        marker="o",
        ls="",
        ms=2,
    )
    ax[1].grid()
    ax[1].set_xlim([Imin, Imax])
    ax[1].set_xlabel(r'$I_{stim} = I_{0}$ ($\mu$A/cm$\,^2$)')
    ax[1].set_ylabel(r'$v$ (mV)')
    ax[1].legend(frameon=False, loc='upper right', bbox_to_anchor = (1.02,1.05), labelspacing=0.1)
    fig.tight_layout()
    fig.savefig('02_hh3Dm.pdf')
    print(np.nanmax(FS))
    plt.show()


if __name__ == "__main__":
    STABLE = 1
    UNSTABLE = 2
    STABLE_LC = 3
    UNSTABLE_LC = 4

    read_xppdat(bifurcation_file, ISI_file)
