import numpy as np
import cPickle as cp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import hera_sim
import hera_cal

plt.style.use(['paper'])
colors = sns.color_palette('colorblind')

with open('lowcadcal_vs_subredcal_ants.cp', 'r') as fp:
     [subsol, subchi, lowsol, lowchi] = cp.load(fp)

hexrange = np.sort(subsol.keys())
Nsim = 256
Nants_list = [] 

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.17, right=0.98, top=0.85, bottom=0.15)

hexNum = 11
#for i,hexNum in enumerate(hexrange):

antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
Nants = len(antpos)
reds = hera_cal.redcal.get_reds(antpos)
Nbls = np.sum([len(s) for s in reds])

NlogN = Nants * np.log2(Nants)

Nants_list.append(Nants)
fraction_tint = np.sort(lowchi[hexNum].keys())
fraction_bls  = np.sort(subchi[hexNum].keys())

chi_sub = []
chi_low = []

for f, tint in zip(fraction_bls, fraction_tint):
    print f, tint
    chi_sub.append(np.mean(subchi[hexNum][f]))
    chi_low.append(np.mean(lowchi[hexNum][tint]))

ax.plot(fraction_bls, chi_sub, 's', ls='--',
        color='k', label='subredcal')
        #color=colors[i], label='subredcal')
ax.plot(fraction_bls, chi_low, '^', ls='-',
        color='k', label='lowcadcal')    
        #color=colors[i], label='lowcadcal')    

ax.axvline(x=1.0*NlogN/Nbls, ls=':', alpha=0.5, color='k')

    #gvar_sub = []
    #gvar_low = []

    #for f,tint in zip(fraction_bls, fraction_tint):
    #    gain_sub = np.zeros([Nants, Nsim], dtype=np.complex)
    #    gain_low = np.zeros([Nants, Nsim], dtype=np.complex)

    #    for ant in range(Nants):
    #        gain_sub[ant] = subsol[hexNum][f][(ant,'Jxx')]
    #        gain_low[ant] = lowsol[hexNum][tint][(ant, 'Jxx')]

    #    gvar_sub.append(np.mean(np.abs(np.diag(np.cov(gain_sub)))))
    #    gvar_low.append(np.mean(np.abs(np.diag(np.cov(gain_low)))))

    #ax.plot(fraction_bls, gvar_sub, 's', ls='--', 
    #        color=colors[i], label='subredcal')
    #ax.plot(fraction_bls, gvar_low, '^', ls='-',
    #        color=colors[i], label='lowcadcal')

#h1, = ax.plot([0],[0], ls = '--', marker='s', color='k', label='subredcal')
#h2, = ax.plot([0],[0], ls = '-',  marker='^', color='k', label='lowcadcal')

ax.legend(fancybox=True)

ax.set_xlim(0.04,0.22)
ax.set_xlabel('Fraction of baselines')

axt = ax.twiny()
axt.set_xlim(ax.get_xlim())
axt.set_xticks([0.05, 0.1, 0.15, 0.20])
ftint = np.sort([float(x) for x in lowchi[hexNum].keys()])
axt.set_xticklabels([0.05, 0.1, 0.15, 0.20])
axt.set_xlabel('Fractional Integration Time', labelpad=15)


ax.set_ylabel(r'$\chi^2_r$')
#ax.set_ylim([1.0,1.6])

ax.grid(which='both', ls='dotted', alpha=0.5)

#Nants_list = np.asarray(Nants_list, dtype=np.float)
#cax = fig.add_axes([0.92,0.15,0.015,0.75])
#
#cmap = mpl.colors.ListedColormap(colors[:len(Nants_list)])
#bounds = np.asarray([28, 94, 172, 274, 400])
##bounds = np.asarray([25, 49, 76, 109, 148, 193, 241])
#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
#                               norm=norm,
#                               ticks=Nants_list,
#                               boundaries= bounds,
#                               spacing='uniform',
#                               orientation='vertical')
#cb.ax.set_title('$N_{\mathrm{ants}}$')

fig.savefig('lowcadcal_vs_subredcal.pdf')

plt.show()
