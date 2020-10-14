import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cPickle as cp
#import _pickle as cp
import pickle as cp
import seaborn as sns
#import hera_sim
import hera_cal
import matplotlib
from matplotlib.patches import Polygon

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

colors = sns.color_palette('colorblind')
plt.style.use(['paper'])

NSIM = 256
steps = np.logspace(0, 6, base=2, num=7, dtype=np.int)

pos = np.loadtxt('antenna_positions_37.dat')
antpos = {k:v for k,v in zip(range(37), pos)}
Nants = len(antpos)
reds = hera_cal.redcal.get_reds(antpos)
Nbls = 1.0*len([pair for subbls in reds for pair in subbls])
sigma = 0.01

w = 0.05
width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)


with open('gain_scaling_bls_legacy_subredcal.cp','r') as fp:
    subsol, num_unibls = cp.load(fp)

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16, right=0.96, top=0.96, bottom=0.15)

num_redbls = []
for num in num_unibls:
    num_redbls.append(np.sum([len(s) for s in reds[0:num]]))
num_redbls = np.take(num_redbls, [0,1,2,3,4,5,6,8,10,13,17,21,25,34,59])
num_unibls = np.take(num_unibls, [0,1,2,3,4,5,6,8,10,13,17,21,25,34,59])

#corner_ant_nums = [0,3,15,21,33,36]
corner_ant_nums = [0, 1, 2, 3, 4, 8, 9, 14, 15, 21, 22, 27, 28, 32, 33, 34, 35, 36]
Nants_corner = len(corner_ant_nums)

# Plot two errorbars- one with all the antennas and 
# another with the corner antennas excluded. 
# Use the same mean-gain-variance for both.
for num in num_redbls:

    # All antennas included
    g = np.zeros([Nants, NSIM], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant]= subsol[num][(ant,'Jxx')]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    avg = np.mean(g_var)
    yerr_all = 0.5*np.std(g_var)

    bplot = ax.boxplot(g_var, positions=[num/Nbls], widths=width([num/Nbls],w),
                       showfliers=False, showmeans=True, meanline=True, 
                       meanprops  =dict(linestyle='-', color='k'),
                       medianprops=dict(linestyle='-', color='w'))
    bplot['medians'][0].set_alpha(0)

    # Without the 6 corner antennas
    g = np.zeros([Nants-Nants_corner, NSIM], dtype=np.complex)
    i = 0
    for ant in range(Nants):
        if ant in corner_ant_nums: continue
        else: 
           g[i]= subsol[num][(ant,'Jxx')]
           i += 1
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    yerr_exc = 0.5*np.std(g_var)

    bplot = ax.boxplot(g_var, positions=[num/Nbls], widths=width([num/Nbls],w),
                       showfliers=False, showmeans=True, meanline=True,
                       meanprops=dict(linestyle='-', color='yellow'))
    box = bplot['boxes'][0]
    box_coords = np.transpose(np.stack([box.get_xdata(), box.get_ydata()]))
    ax.add_patch(Polygon(box_coords, facecolor=colors[-3]))
    bplot['medians'][0].set_alpha(0)

    print '{0:3d}\t{1:.4f}\t{2:.4f}\t{3:2.2f}'.format(num, yerr_all, 
          yerr_exc, (yerr_all - yerr_exc)/yerr_all * 100)

    #ax.errorbar(num/Nbls, y= avg, yerr= yerr_exc,
    #            fmt='none', color='k', ms=3,
    #            elinewidth=7, capsize=0, 
    #            ecolor=colors[-3])
    #
    #ax.errorbar(num/Nbls, y= avg, yerr= yerr_all,
    #            fmt='none', color= 'k', ms=3)


num_redbls = np.asarray(num_redbls, dtype=np.float)
num_unibls = np.asarray(num_unibls, dtype=np.float)

ax.plot(num_redbls/Nbls, 0.55*(Nants/num_redbls),'-', 
        color=colors[0])
ax.plot(num_redbls[:11]/Nbls, 2.6/num_unibls[:11]**1.5, '--', 
        color=colors[2])


ax.annotate(r'$\propto N^{-1}_{\rm{bl;ant}}$', 
            xy=(0.7,0.08),  color='k',
            fontweight=990, fontsize=18)

ax.annotate(r'$\propto N^{-1.5}_{\rm{ubl;ant}}$', 
            xy=(0.25,0.4), color='k',  
            fontweight=990, fontsize=18)

ax.loglog()
ax.set_xlabel(r'Fraction of all baselines')
ax.set_ylabel(r'$\sigma_{g}^2 $')
ax.grid(which='both', ls='dotted', alpha=0.7)
ax.set_xlim([0.09, 1.3])
ax.set_xticks([0.1,0.5,1])
ax.set_xticklabels([0.1,0.5,1])

#fig.savefig('subredcal_gain_scaling_with_bls.pdf')
#
#
## --------------------------
## For different array sizes
## --------------------------
#with open('gain_scaling_bls_ants_subredcal_omnical_LN.cp','r') as fp:
#    subsol = cp.load(fp)
#
#Nsim = 256
#sigma = 0.1
#
#Nants_list = np.sort(subsol.keys())
#
#fig,ax = plt.subplots(1,1)
#fig.subplots_adjust(left=0.16, right=0.85, top=0.99, bottom=0.14)
#
#w = 0.05
#width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
#
#for i,Nants in enumerate(Nants_list):
#
#    Nbls = subsol[Nants]['Nbls']
#    NlogN = subsol[Nants]['NlogN']
#    f_NlogN = 2.0*NlogN/Nbls
#
#    fraction_bls = filter(lambda x: type(x)!=str, subsol[Nants].keys())
#    fraction_bls = np.sort(fraction_bls)
#
#    num_redbls = subsol[Nants]['num_redbls']
#    num_unibls = subsol[Nants]['num_unibls']
#    y1 = []; y2 = [];
#
#    g = np.zeros([Nants, Nsim], dtype=np.complex)
#
#    for f in fraction_bls:
#        for ant in range(Nants):
#            g[ant] = subsol[Nants][f][(ant,'Jxx')]
#
#        g_var = np.abs(np.diag(np.cov(g))) / sigma**2
#
#        bplot = ax.boxplot(g_var, positions=[f], patch_artist=True, 
#                           widths=width([f],w), showfliers=False,
#                           medianprops=dict(color='w'), showmeans=True, 
#                           meanline=True,
#                           meanprops=dict(linestyle='-', color='k'))
#        bplot['boxes'][0].set_facecolor(colors[i])
#        bplot['medians'][0].set_alpha(0)
#        
#
#        #avg = np.mean(g_var)
#        #y1.append(avg+ 0.5*np.std(g_var))
#        #y2.append(avg- 0.5*np.std(g_var))
#        #ax.plot(f, avg, 'o', color=colors[i])
#
#    #ax.plot(f, avg, 'o', color=colors[i])
#    #ax.fill_between(fraction_bls, y1=y1, y2=y2, color=colors[i], 
#    #                alpha=0.3, label='%d'%Nants)
#
#    num_redbls = np.sort(np.asarray(num_redbls, dtype=np.float))
#    num_unibls = np.sort(np.asarray(num_unibls, dtype=np.float))
#
#    print (Nants)
#    print (num_redbls)
#    print (num_unibls)
#
#    ax.plot(num_redbls[:5]/Nbls, 2.5/num_unibls[:5]**1.5, '--',
#            color=colors[i], alpha=0.7)
#    ax.plot(f_NlogN, 0.52*Nants/(2*NlogN), 'X', 
#            color='k', markersize=7)
#    ax.plot(num_redbls/Nbls, 0.52*Nants/num_redbls, '-',
#            color=colors[i], alpha=0.7)
#
#ax.loglog()
#ax.grid(which='both', ls='dotted', alpha=0.7)
#Nants_list = np.asarray(Nants_list, dtype=np.float)
#cax = fig.add_axes([0.90,0.15,0.015,0.75])
#
#cmap = mpl.colors.ListedColormap(colors[:len(Nants_list)])
#bounds = np.asarray([160, 274, 481, 886, 1400])
#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
#                               norm=norm,
#                               ticks=Nants_list,
#                               boundaries= bounds,
#                               spacing='uniform',
#                               orientation='vertical')
#cb.ax.set_title('$N_{\mathrm{ants}}$')
#
#ax.loglog()
#ax.set_xlabel('Fraction of all baselines', fontsize=20)
#ax.set_ylabel(r'$\sigma_{g}^2 $', fontsize=20)
#ax.grid(which='both', ls='dotted', alpha=0.5)
#ax.set_xlim([3e-3,1.2])
#ax.set_xticks([0.01, 0.1, 1])
#ax.set_xticklabels([0.01, 0.1, 1])
#
#fig.savefig('subredcal_gain_scaling_bls_ants.pdf')

plt.show()
