import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp
import seaborn as sns

plt.style.use(['paper'])
colors = sns.color_palette('colorblind')

Nsim = 4096
sigma = 0.01
sigma_low = 0.05
steps = np.logspace(0, 10, base=2, num=11, dtype=np.int)
fbls = 0.04

with open('nbls_per_vis_hexNum4-20.cp','r') as fp:
    Nbls = cp.load(fp)

Nbls_hexNum9 = np.sum([v for v in Nbls[9].values()])

with open('lowcadcal_gains_mvis_converge_hexNum9.cp','r') as fp:
    lowsol, mvis_low, tg_low, tv_low, Nants = cp.load(fp)

with open('subredcal_gains_mvis_converge.cp','r') as fp:
    subsol, mvis_sub, tg_sub, tv_sub, Nants = cp.load(fp)


## Gain convergence for both subredcal and lowcadcal

#fig,ax = plt.subplots(2,1, figsize=(6.4, 6.4))
#fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.15)
fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.15)

low_g_matr = np.zeros([Nants, Nsim], dtype=np.complex)
tru_g_mlow = np.zeros([Nants, Nsim], dtype=np.complex)
sub_g_matr = np.zeros([Nants, Nsim], dtype=np.complex)
tru_g_msub = np.zeros([Nants, Nsim], dtype=np.complex)

for ant in range(Nants):
    low_g_matr[ant] = lowsol[(ant,'Jxx')]
    tru_g_mlow[ant] = tg_low[(ant,'Jxx')]
    sub_g_matr[ant] = subsol[(ant,'Jxx')]
    tru_g_msub[ant] = tg_sub[(ant,'Jxx')]

for nsim in steps:

    low_gavg = np.mean(low_g_matr.reshape(Nants,-1,nsim), axis=2)
    tru_glow = np.mean(tru_g_mlow.reshape(Nants,-1,nsim), axis=2)

    dev = np.abs(low_gavg - tru_glow)**2 / sigma_low**2 * Nants

    h1 = ax.errorbar(nsim, y=np.mean(dev), yerr= 0.5*np.std(dev),
                     fmt='^', color=colors[2], label='lowcadcal',
                     fillstyle='none', markersize=7, mew=2, alpha=0.8)

    sub_gavg = np.mean(sub_g_matr.reshape(Nants,-1,nsim), axis=2)
    tru_gsub = np.mean(tru_g_msub.reshape(Nants,-1,nsim), axis=2)

    dev = np.abs(sub_gavg - tru_gsub)**2 / sigma**2 * (fbls * Nbls_hexNum9 / Nants)

    h2 = ax.errorbar(nsim, y=np.mean(dev), yerr= 0.5*np.std(dev),
                     fmt='s', color=colors[1], label='subredcal',
                     fillstyle='none', markersize=7, mew=2, alpha=0.8)

ax.plot(steps, 1.0/(steps), '--k', alpha=0.7)

ax.loglog()
ax.grid(ls='dotted', which='both', alpha=0.3)
ax.set_xlabel('Independent noise realizations averaged over', fontsize=18)
ax.set_ylabel(r'$\left|\bar{{\bf g}}-{\bf g}^{\mathrm{input}}\right|^2/\sigma_g^2$', fontsize=18)
ax.legend(handles=[h1,h2], fancybox=True)


### Visibility convergence for subredcal and lowcadcal
#
#k = (0,11,'xx')
#
#for nsim in steps:
#    low_vis = np.mean(np.reshape(mvis_low[k], (-1,nsim)), axis=1)
#    low_tru = np.mean(np.reshape(tv_low[k],   (-1,nsim)), axis=1)
#
#    sub_vis = np.mean(np.reshape(mvis_sub[k], (-1,nsim)), axis=1)
#    sub_tru = np.mean(np.reshape(tv_sub[k],   (-1,nsim)), axis=1)
#    
#    dev = np.abs(low_vis-low_tru)**2/sigma**2 * Nbls[9][k]
#    
#    h1, = ax[1].plot(nsim, np.mean(dev), #yerr= 0.5*np.std(dev), 
#                     marker='^', color=colors[2], label= 'lowcadcal',
#                     fillstyle='none', markersize=7, mew=2, alpha=0.8)
#
#    dev = np.abs(sub_vis-sub_tru)**2/sigma**2 * Nbls[9][k]
#    
#    h2, = ax[1].plot(nsim, np.mean(dev), #yerr= 0.5*np.std(dev), 
#                     marker='s', color=colors[1], label= 'subredcal',
#                     fillstyle='none', markersize=7, mew=2, alpha=0.8)
#
#ax[1].plot(steps, 1.0/(steps), '--k', alpha=0.7)
#    
#ax[1].loglog()
#ax[1].grid(ls='dotted', which='both', alpha=0.3)
#ax[1].set_xlabel('Independent noise realizations averaged over', fontsize=16)
#ax[1].set_ylabel(r'$\left|\bar{V_{\alpha}}-V_{\alpha}^{\mathrm{input}}\right|^2/\sigma_V^2$', fontsize=18)
##ax[1].legend(handles=[h1,h2], fancybox=True)

fig.savefig('gain_vis_convergence.pdf')

plt.show()
