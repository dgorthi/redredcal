import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import cPickle as cp
from collections import OrderedDict
import matplotlib as mpl

NSIM = 256
steps = np.logspace(0, 6, base=2, num=7, dtype=np.int)

plt.style.use(['paper'])
#plt.style.use(['presentation'])
line_style=[':','--','-']
colors = sns.color_palette('colorblind')

ants = np.loadtxt('antenna_positions_37.dat')
antpos = {k:v for k,v in zip(range(len(ants)), ants)}
Nants = len(antpos)

w = 0.05
width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)

with open('gain_scaling_snr_lowcadcal.cp','r') as fp:
    lowsol = cp.load(fp)

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16, right=0.95, top=0.98, bottom=0.16)

for sgl,sol in lowsol.items():
    g = np.zeros([len(antpos), NSIM], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant]= sol[(ant,'Jxx')]

    g_var = np.abs(np.diag(np.cov(g)))
    avg = np.mean(g_var)

    bplot = ax.boxplot(g_var, positions=[1.0/sgl], widths=width([1.0/sgl],w),
                       showfliers=False, showmeans=True, meanline=True,
                       meanprops  =dict(linestyle='-', color='k'),
                       medianprops=dict(linestyle='-', color='w'))
    bplot['medians'][0].set_alpha(0)

    #ax.errorbar(1.0/sgl, y= avg, yerr= np.std(g_var), 
    #            fmt='none', color='k', ms=3, capsize=5)

snr = np.linspace(0.2, 5, num=200)
ax.plot(snr, 1.1/(Nants*snr**2),'--',color='k')
ax.annotate(r'$\propto$ $\frac{1}{\rm{SNR}^2}$', xy=(2, 0.02),
            fontsize=24)

ax.loglog()
ax.set_xlabel(r'\rm{SNR of} $V_{ij}^{M}$')
ax.set_ylabel(r'$\sigma_{g}^2 $')
ax.grid(which='both', ls='dotted', alpha=0.7)
ax.set_xlim([0.2,5])
ax.set_xticks([0.3, 1, 3])
ax.set_xticklabels([0.3, 1, 3])

fig.savefig('lowcadcal_gain_scaling_with_snr.pdf')

## ---------------------------------
##  Function of number of antennas
## ---------------------------------
#
#Nsim = 256
#with open('gain_scaling_snr_ants_lowcadcal_omnical_LN.cp','r') as fp:
##with open('gain_scaling_snr_ants_lowcadcal_lincal.cp','r') as fp:
#    lowsol = cp.load(fp)
#
#Nants = np.sort(lowsol.keys())
##Nants = Nants.take([0,1,3,5])
#sigma_g = {}
#for na in Nants:
#    sigma_g[na] = {}
#
#fig, ax = plt.subplots(1,1)
#fig.subplots_adjust(left=0.16, right=0.95, top=0.98, bottom=0.16)
#
#for i, na in enumerate(Nants):
#    #lowsol[na].pop(10.0)
#    #lowsol[na].pop(6.309573444801933)
#    y1 = []; y2 = [];
#    
#    sigma_range = np.sort(lowsol[na].keys())
#    bplot = []
#
#    for sigma in sigma_range:
#        g = np.zeros([na, Nsim], dtype=np.complex)
#        for ant in range(na):
#            g[ant] = lowsol[na][sigma][(ant,'Jxx')]
#        g_var = np.abs(np.diag(np.cov(g)))
#
#        bplot = ax.boxplot(g_var, positions=[1.0/sigma], patch_artist=True, 
#                           widths=width([1.0/sigma],w), showfliers=False,
#                           medianprops=dict(color=colors[i], linestyle='-'))
#        #bplot['boxes'][0].set_facecolor(colors[i])
#        plt.setp(bplot['boxes'],    color=colors[i]) 
#        plt.setp(bplot['whiskers'], color=colors[i]) 
#        plt.setp(bplot['caps'],     color=colors[i])
#        
#
#        #avg = np.mean(g_var)
#        #y1.append(avg + np.std(g_var))
#        #y2.append(avg - np.std(g_var))
#
#        #ax.errorbar(1.0/sigma, y=avg, yerr = np.std(g_var),
#        #            fmt='none', color=colors[i])
#    #ax.fill_between(1.0/sigma_range, y1=y1, y2=y2, color=colors[i], alpha=0.2)
#
#Nants = np.asarray(Nants, dtype=np.float)
#
#for i, na in enumerate(Nants):
#    ax.plot(1.0/sigma_range, (na**-1)*(sigma_range**2),
#            ls = '-', color=colors[i], label=int(na))
#
#ax.annotate(r'$\propto$ $\frac{1}{(\rm{SNR})^2\; N_{\rm{ants}}}$', xy=(1.9, 0.002),
#            fontsize=24)
#
##cax = fig.add_axes([0.90,0.15,0.015,0.75]) 
##cmap = mpl.colors.ListedColormap(colors[:len(Nants)])
##bounds = np.asarray([160, 274, 481, 886, 1400])
##norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
##cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
##                               norm=norm,
##                               ticks=Nants,
##                               boundaries= bounds,
##                               spacing='uniform',
##                               orientation='vertical')
##cb.ax.set_title(r'$N_{\rm{ants}}$')
#
#ax.loglog()
#ax.set_xlabel(r'SNR of $V_{ij}^{M}$', fontsize=18)
#ax.set_ylabel(r'$\sigma_{g}^2 $', fontsize=18)
#ax.grid(which='both', ls='dotted', alpha=0.5)
#ax.set_xlim([0.2,12])
#ax.set_xticks([0.3,1,10])
#ax.set_xticklabels([0.3,1,10])
#ax.legend(fancybox=True, loc=3)
#
#fig.savefig('lowcadcal_gain_scaling_snr_ants.pdf')

plt.show()
