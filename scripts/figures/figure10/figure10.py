import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors

def map_plot(Amin, Amax, fMin, fMax, I_0, Alen, flen, model):

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 18
    plt.rcParams["mathtext.fontset"] = "stix"

    # loading data
    fname = "average_map_"+str(model)+"_Amin_"+str(Amin)+"_Amax_"+str(Amax)+"_fmin_"+str(fMin)+"_fmax_"+str(fMax)+"kHz_I0_"+str(I_0)+"_"+str(Alen)+"_"+str(flen)+""
    file = ""+fname+".csv"
    raw_data = np.loadtxt(file, delimiter=',')

    stim_current = raw_data[0,1:]
    stim_freq = raw_data[1:, 0]
    delta_v = raw_data[1:,1:]

    del raw_data

    # construct colormap
    N = 200
    vmin = np.amin(delta_v)
    vminabs = np.abs(vmin)
    vmax = np.amax(delta_v)

    print()
    bottom = plt.cm.Blues(np.linspace(1,0.1, int(N*(vminabs/(vminabs+vmax)))))
    top = plt.cm.Oranges(np.linspace(0.1, 1, int(N*(vmax/(vminabs+vmax)))))


    newcolors = np.vstack((bottom, top))
    newcmp = mcolors.LinearSegmentedColormap.from_list('my_colormap', newcolors)

    A_stim = np.array(stim_current)
    f_stim_4 = A_stim / (2 * np.pi * 23) + 0.02 # f_stim = A_stim / (2π * 50)
    # plot
    fig, ax = plt.subplots(figsize=(6,5))
    map = plt.pcolormesh(stim_current, stim_freq, delta_v, cmap=newcmp, rasterized=True)
    
   # plt.plot(A_stim, f_stim_4, linestyle='--', linewidth = 2, color='black')
    # Ajustar o range do gráfico (a) para manter o intervalo original
    ax.set_xlim(stim_current.min(), stim_current.max())
    ax.set_ylim(stim_freq.min(), stim_freq.max())

    #ax.yaxis.set_inverted(True)
    # plt.contour(stim_current, stim_freq, delta_v, 0, linewidths=0.5, colors=['k'])
    cbar = plt.colorbar(map)
    cbar.set_label("$\\langle \\Delta v \\rangle$")
    plt.ylabel(r'stimulation frequency, $\: f_{\mathrm{stim}}$ (kHz)')
    plt.xlabel(r'stimulation amplitude, $\: A_{\mathrm{stim}}$ ($\mathrm{\mu A / cm^2}$)')
    #plt.title(f"Constant Stimulation: $I_0 = {{{int(I_0)}}}\: \mu A / cm^2$", fontsize=16)
    plt.tight_layout()
    plt.savefig(fname+'.pdf', dpi=400)
    plt.show()

models = [
    ["WangBuzsaki", 10, 1000, 100, 0.1, 10.0, 100, 10, 200],
    ["TraubMiles", 10, 3000, 100, 0.1, 10.0, 100, 10, 200],
    ["Erisir", 10, 1400, 100, 0.1, 5.0, 100, 10, 200],
    ["TermanGPe", 10, 600, 100, 0.1, 3.0, 100, 10, 500],
    ["TermanSTN", 10, 2000, 100, 0.1, 10.0, 100, 10, 500],
    ["HH", 10, 7000, 100, 0.1, 3.5, 100, 0, 425],
    ["Rgc1", 10, 4000, 100, 0.5, 20.0, 100, 10, 100],
    ["Rgc2", 10, 4000, 100, 0.5, 20.0, 100, 10, 100]
    ]
for modelname, Amin, Amax, Alen, fmin, fmax, flen, I0, Alim in models:
    map_plot(Amin, Amax, fmin, fmax, I0, Alen, flen, modelname)