import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 24
plt.rcParams["mathtext.fontset"] = "stix"

params = [
      ("ErisirN2", 10, 1400, 100, 5000, 10, 300, 300, 1e-2, 0.22),
      ("TraubMiles", 10, 3000, 100, 7500, 10, 300, 300, 1e-2,  1),
      ("WangBuzsaki", 10, 500,100,5000,10,500, 500, 1e-2, 0.5),
      ("TermanGP", 10, 1000, 800, 5000, 10, 400, 400,1e-2, 0.8),
      ("TermanSTN", 10, 1500, 100, 10000, 10, 400, 400,1e-2, 0.45),
      ("rgcTypeI", 10, 4000, 500, 20000, 10, 300, 300, 1e-2, 0.5),
      ("rgcTypeII", 10, 4000, 100, 20000, 10, 300, 300, 1e-2, 0.5),
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
    fminplot = params[i][9]


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

    fig, axs = plt.subplots(figsize=(5,4))
    fname = 'dynamics-'+str(Model)+'/periods_Amin_'+str(Amin)+'_Amax_'+str(Amax)+'_fMin_'+str(fmin)+'_fMax_'+str(fmax)+'_I0_'+str(I0)+'_size_'+str(lengthA)+'_x_'+str(lengthf)+'.csv'
    raw_data = np.loadtxt(fname, delimiter='\t')
    stim_current = raw_data[0,1:]
    stim_freq = raw_data[1:, 0]
    periods = raw_data[1:,1:]
    periods = np.ma.masked_equal(periods, -2.0)
    psm = axs.pcolormesh(stim_current,stim_freq, periods, cmap=cmap, rasterized=True, vmin=-2, vmax=7)
   # axs.set_xlim([10,1000])
    plt.ylim([fminplot,fmax/1000])


    # Calcular as funções para as linhas tracejadas
    A_stim = np.array(stim_current)  # Eixo x
    
    f_stim_25 = A_stim / (2 * np.pi * 23) + 0.02 # f_stim = A_stim / (2π * 50) STN
    f_stim_15 = A_stim / (2 * np.pi * 15.1)  # f_stim = A_stim / (2π * 50) WB
   

    plt.savefig('periods-'+str(Model)+'_Amin_'+str(Amin)+'_Amax_'+str(Amax)+'_fMin_'+str(fmin)+'_fMax_'+str(fmax)+'_I0_'+str(I0)+'_size_'+str(lengthA)+'_x_'+str(lengthf)+'.pdf', dpi=300)
    plt.show()

