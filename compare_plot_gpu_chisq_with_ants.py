import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['lines.markersize'] = 6
plt.style.use(['paper'])
colors = sns.color_palette('colorblind')

#with open('scale_noise_wchisq_hexNum4-23.cp','r') as fp:
#    gains, model_vis, chisq, true_gains, true_vis, Nants_list = cp.load(fp)

with open('compare_gainvar_chisq_hexNums4-40_final.pkl','rb') as fp:
    gains, chisq, Nubls = pickle.load(fp)

Nants_list = gains.keys()

fig, ax = plt.subplots(2,1, gridspec_kw={'height_ratios': [3, 1]}, 
                      figsize=(6.4,5.8), sharex=True)
fig.subplots_adjust(left=0.15,right=0.98,bottom=0.15,top=0.95)

for nants in Nants_list:
    h3 = ax[0].errorbar(nants, y = chisq[nants]['redcal']['avg'], 
                        yerr=  0.5*chisq[nants]['redcal']['std'], 
                        fmt='o', color='k', label='redcal',
                        ms=5, capsize=5, fillstyle='none')

    h1 = ax[0].errorbar(nants, y = chisq[nants]['subredcal']['avg'], 
                        yerr=  0.5*chisq[nants]['subredcal']['std'],
                        fmt='s', color=colors[1], label='subredcal',
                        ms=5, capsize=5)

    h2 = ax[0].errorbar(nants, y = chisq[nants]['lowcadcal']['avg'], 
                        yerr=  0.5*chisq[nants]['lowcadcal']['std'], 
                        marker='^', color=colors[2], label='lowcadcal',
                        ms=5, capsize=5)
                   
    ax[0].errorbar(nants, y = chisq[nants]['subredcal_sub']['avg'], 
                   yerr=  0.5*chisq[nants]['subredcal_sub']['std'], 
                   fmt='s', color=colors[1], alpha=0.6, 
                   ms=5, capsize=5, fillstyle='none')

    ax[0].errorbar(nants, y = chisq[nants]['lowcadcal_low']['avg'], 
                   yerr=  0.5*chisq[nants]['lowcadcal_low']['std'], 
                   fmt='^', color=colors[2], alpha=0.6,
                   ms=5, capsize=5, fillstyle='none')

    #redbls = [bl for bl in gains[hexNum]['subredcal'].keys() if bl[-1]=='xx']
    #Nvar.append(len(redbls))


# Get the number of unique baselines in the subredcal calibrator
# This is a proxy for the calibrator size
#import hera_sim
#import hera_cal
#
#Nubls_subredcal = []
#Nants_list = []
#
#hexNums = list(np.arange(4,25)) + [25,27,30,33,35,37,40]
#
#for hexNum in hexNums:
#    antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
#    reds   = hera_cal.redcal.get_reds(antpos)
#
#    Nants = len(antpos)
#    Nbls  = np.sum([len(s) for s in reds]) 
#    Nants_list.append(Nants)
#
#    NlogN = int(Nants * np.log2(Nants))
#    N = 3
#    subreds = reds[0:N]
#    Ncorr = np.sum([len(s) for s in subreds])
#    while (Ncorr < NlogN):
#        N = N+1
#        subreds = reds[0:N]
#        Ncorr = np.sum([len(s) for s in subreds])
#
#    N = N-1
#    subreds = reds[0:N]
#    Nubls_subredcal.append(len(subreds))

ax[0].semilogx()
ax[0].set_xlim([50, 1.5e4])
ax[0].set_ylim([0.9,1.4])
ax[0].set_xticklabels([])
ax[0].set_ylabel(r'$\chi_r^2$')
ax[0].legend(handles=[h1,h2,h3], handlelength=0, fancybox=True, ncol=3)
ax[0].grid(ls='dotted', which='both', alpha=0.6)

Nvar = [v for v in Nubls.values()]
Nants_list = [k for k in Nubls.keys()]
ax[1].semilogx(Nants_list[1:], Nvar[1:], '--o', color='k')

ax[1].set_xlim([50, 1.5e4])
ax[1].set_xlabel('Number of antenna elements')
ax[1].set_ylim([5,15])
ax[1].set_yticks([8,11,14])
ax[1].set_ylabel('$N_{\mathrm{ubl}}$')
ax[1].grid(ls='dotted', which='both', alpha=0.6)
fig.savefig('chi_scaling_antennas.pdf')

plt.show()  

