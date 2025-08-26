import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams["mathtext.fontset"] = "stix"

params = [
      # ("ErisirN2", 10, 1400, 100, 5000, 10, 300, 300, 1e-2, 0.22),
      # ("TraubMiles", 10, 3000, 100, 7500, 10, 300, 300, 1e-2,  1),
      #  ("WangBuzsaki", 10, 500,100,5000,10, 500, 500, 1e-2, 0.5),
     ("HodgkinHuxley", 10, 7000, 400, 3500, 0, 900, 900, 1e-2, 0.4),
      # ("TermanGP", 500, 750, 4000, 4500, 10, 400, 400,1e-2, 4.0),
      # ("TermanSTN", 10, 1500, 100, 10000, 10, 400, 400,1e-2, 0.45),
      # ("rgcTypeI", 10, 4000, 500, 20000, 10, 300, 300, 1e-2, 0.5),
      # ("rgcTypeII", 10, 4000, 100, 20000, 10, 300, 300, 1e-2, 0.5)
]

for i in range(len(params)):
    Model = params[i][0]
    Amin = params[i][1]
    Amax = params[i][2]
    fmin = params[i][3]
    fmax = params[i][4]
    I0 = params[i][5]
    lengthA = params[i][6]
    lengthf = params[i][7]
    fname = 'lyapunov_Amin_'+str(Amin)+'_Amax_'+str(Amax)+'_fMin_'+str(fmin)+'_fMax_'+str(fmax)+'_I0_'+str(I0)+'_size_'+str(lengthA)+'_x_'+str(lengthf)+'.csv'
    raw_data = np.loadtxt(fname, delimiter='\t')

    stim_current = raw_data[0,1:]
    stim_freq = raw_data[1:, 0]
    lyapunov_exponents = raw_data[1:,1:]

    # construct colormap
    N = 500
    vmin = np.amin(lyapunov_exponents)
    vminabs = np.abs(vmin)
    vmax = np.amax(lyapunov_exponents)


    colorsneg = [(232/255,239/255,249/255), (0,0,90/255)]
    colormapNegative = mcolors.LinearSegmentedColormap.from_list('negative_colors', colorsneg)
    bottom = colormapNegative(np.linspace(0,1, int(N*(vminabs/(vminabs+vmax)))))
    colorspos = [(248/255,244/255,180/255),(139/255,2/255,2/255)]
    colormapPositive = mcolors.LinearSegmentedColormap.from_list('positive_colors', colorspos, N = int(N*(vmax/(vminabs+vmax))))
    top = colormapPositive(np.linspace(0,1, int(N*(vmax/(vminabs+vmax)))))

    colorszero = [(0.0,0.0,0.0), (0.0,0.0,0.0)]
    colormapZero = mcolors.LinearSegmentedColormap.from_list('zero', colorszero)
    zeros = colormapZero(np.linspace(0,1,2))

    epsilon = params[i][8]

    lpos = np.ma.masked_less(lyapunov_exponents, epsilon)
    lzero = np.ma.masked_outside(np.abs(lyapunov_exponents), 10**(-20), epsilon)
    lneg = np.ma.masked_greater(lyapunov_exponents, -epsilon)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')


    pos = axs[0].pcolormesh(stim_current, stim_freq, lpos, cmap=colormapPositive, rasterized = True)
    axs[0].pcolormesh(stim_current, stim_freq, lzero, cmap= colormapZero, rasterized = True)
    axs[0].pcolormesh(stim_current, stim_freq, lneg, cmap=colormapNegative, rasterized = True)
    axs[1].tick_params(labelleft=False)

    newcolors = np.vstack((bottom, zeros, top))
    newcmp = mcolors.LinearSegmentedColormap.from_list('my_colormap', newcolors)
    norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax)

    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=newcmp),
                ax=axs[0], location = 'left')

    cbar.set_label('Maximum $\lambda$')
    #cbar.ax.set_title('Maximum $\lambda$')
    axs[0].set_ylabel('$f_{\mathrm{stim}}$ (kHz)')
    #axs[0].set_xlabel('$A_{stim}$ ($\mu$A/cm²)')
    fig.supxlabel('$A_{\mathrm{stim}}$ ($\mu$A/cm²)')


    cmap = ListedColormap([
        (0.85, 0.85, 0.85),  # Não-spike (cinza prateado)
        (1.0, 0.6, 0.4),     # Caótico (laranja avermelhado vívido)
        (0.1, 0.1, 0.1),     # Quasi-periódico (grafite escuro)
        (0.6, 0.8, 1.0),     # Período 1 (azul pastel mais escuro)
        (0.3, 0.75, 0.4),    # Período 2 (verde folha vibrante)
        (0.5, 0.3, 0.8),     # Período 3 (roxo mais saturado)
        (0.85, 0.3, 0.6),    # Período 4 (rosa magenta vivo)
        (1.0, 0.85, 0.2),    # Período 5 (amarelo ouro vivo)
        (0.4, 0.0, 0.4),    # Período 6 (vinho profundo com mais roxo)
        (0.2, 0.25, 0.8)     # Período ≥7 (azul algo mais sóbrio)
    ])





    fname = 'periods_Amin_'+str(Amin)+'_Amax_'+str(Amax)+'_fMin_'+str(fmin)+'_fMax_'+str(fmax)+'_I0_'+str(I0)+'_size_'+str(lengthA)+'_x_'+str(lengthf)+'.csv'
    raw_data = np.loadtxt(fname, delimiter='\t')
    stim_current = raw_data[0,1:]
    stim_freq = raw_data[1:, 0]
    periods = raw_data[1:,1:]
    periods = np.ma.masked_equal(periods, -2.0)
    psm = axs[1].pcolormesh(stim_current,stim_freq, periods, cmap=cmap, rasterized=True, vmin=-2, vmax=7)
    cbar = fig.colorbar(psm, ax=axs[1])
    cbar.set_ticks([-2, -1,  0, 1, 2, 3, 4, 5, 6, 7])
    cbar.set_ticklabels([r"NS", r"Chaotic", r"QP", r"$P_1$",r"$P_2$",r"$P_3$",r"$P_4$",r"$P_5$",r"$P_6$",r"$P_{\geq 7}$"])
    cbar.set_label('Dynamics')

    for label in cbar.ax.get_yticklabels():
        label.set_verticalalignment('center')

    #Adicionar rótulos (a) e (b) dentro do triângulo inferior direito de cada figura
    axs[0].text(0.95, 0.02, '(a)', transform=axs[0].transAxes, fontsize=16, 
                va='bottom', ha='right', color='black')
    axs[1].text(0.95, 0.02, '(b)', transform=axs[1].transAxes, fontsize=16,
                va='bottom', ha='right', color='black')

    # Calcular as funções para as linhas tracejadas
    A_stim = np.array(stim_current)  # Eixo x
    f_stim_50 = A_stim / (2 * np.pi * 40)  # f_stim = A_stim / (2π * 50)
    f_stim_425 = A_stim / (2 * np.pi * 425)  # f_stim = A_stim / (2π * 425)

    # Adicionar as linhas tracejadas no gráfico (a)
    axs[0].plot(A_stim, f_stim_50, linestyle='--', linewidth = 1, color='red', label=r'$A = 40 \, \mathrm{mV}$')
    axs[0].plot(A_stim, f_stim_425, linestyle='--', linewidth = 1, color='orange', label=r'$A = 425 \, \mathrm{mV}$')

    # Ajustar o range do gráfico (a) para manter o intervalo original
    axs[0].set_xlim(stim_current.min(), stim_current.max())
    axs[0].set_ylim(stim_freq.min(), stim_freq.max())

    # Adicionar legendas no triângulo inferior direito (região branca)
    axs[0].legend(
        loc='lower right', fontsize=14, frameon=False,  # Sem fundo ou borda
        bbox_to_anchor=(1.0, 0.03)  # Posicionar manualmente dentro da área branca
    )

    plt.savefig('figure6.pdf', dpi=300)
    plt.show()