import numpy as np
import hera_cal
import hera_sim
import cPickle as cp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

colors = sns.color_palette('colorblind')
plt.style.use(['paper'])

def get_subreds(num_bls, reds):
    """ Get the list of redundant sets that have 
        number of baseslines <= num_bls
    """
    N=0; i=0;
    for num, subbls in enumerate(reds):
        if np.isclose(N, num_bls) or N > num_bls:
           return num
        N += len(subbls)

antpos = hera_sim.antpos.hex_array(4, split_core=False, outriggers=0)
Nants = len(antpos)
reds = hera_cal.redcal.get_reds(antpos)
Nbls = 1.0*np.sum([len(s) for s in reds])

Nsim = 256
sigma = 0.01

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16, right=0.96, top=0.96, bottom=0.15)

with open('gain_scaling_bls_subredcal.cp','r') as fp:
    subsol, uni_bl_types, subvisbls, numbls, ants_per_unibl =cp.load(fp)

total_baselines  = []
unique_baselines = []
uni_bl_per_ant   = []
gain_variance    = []

for num, bltype in enumerate(uni_bl_types):
    if num<3: continue

    total_baselines.append(numbls[bltype])
    unique_baselines.append(num)
    nants = 1.0*np.sum([v for k,v in ants_per_unibl.items() if k in uni_bl_types[:num]])
    uni_bl_per_ant.append(nants/Nants)
    
    g = np.zeros([Nants, Nsim], dtype=np.complex)
    for ant in antpos.keys():
        g[ant] = subsol[bltype][(ant,'Jxx')]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2

    gain_variance.append(np.mean(g_var)) 

total_baselines  = np.asarray(total_baselines, dtype=np.float)
unique_baselines = np.asarray(unique_baselines, dtype=np.float) 
uni_bl_per_ant   = np.asarray(uni_bl_per_ant, dtype=np.float)
gain_variance    = np.asarray(gain_variance)

A = np.transpose(np.vstack((Nants*total_baselines**-1, unique_baselines**-1.5)))
x, res, rank, s = np.linalg.lstsq(A, gain_variance, rcond=None)

uni_bl_types_plot = [0,1,2,3,4,5,6,7,8,9,11,13,16,20,24,28,37,62]

for num, bltype in enumerate(uni_bl_types_plot):
    if num<3: continue

    g = np.zeros([Nants, Nsim], dtype=np.complex)
    for ant in antpos.keys():
        g[ant] = subsol[bltype][(ant,'Jxx')]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2

    ax.errorbar(numbls[bltype]/Nbls, y = np.mean(g_var),
                yerr= 0.5*np.std(g_var), fmt='o', 
                color='k')

ax.plot(total_baselines/Nbls, (x[0]*Nants*total_baselines**-1)+(x[1]*uni_bl_per_ant**-1.5), 
        '-', color=colors[2], lw=2)
ax.plot(total_baselines[3:]/Nbls, 0.53*(Nants/total_baselines[3:]),'-.', 
        color=colors[9], lw=2)
ax.plot(total_baselines[:15]/Nbls, 2.6/uni_bl_per_ant[:15]**1.5, '--', 
        color=colors[6], lw=2)

ax.annotate(r'$12.062 \; N^{-1}_{bls} + 1.902 \; N^{-1.5}_{uv}$', 
            xy=(0.25,0.3),  color='k',
            fontsize=18)

    
ax.loglog()
ax.set_xlabel(r'Fraction of all baselines')
ax.set_ylabel(r'$\sigma_{\left| g\right|}^2 $')
ax.grid(which='both', ls='dotted', alpha=0.7)
ax.set_xlim([0.09, 1.3])
ax.set_xticks([0.1,0.5,1])
ax.set_xticklabels([0.1,0.5,1])

fig.savefig('subredcal_gain_scaling_with_bls_v2.pdf')

plt.show()
