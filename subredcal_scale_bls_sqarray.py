import numpy as np
import hera_sim
import hera_cal
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as cp

colors = sns.color_palette('colorblind')
plt.style.use(['paper'])

def square_array(side, sep=14.6):
    """
    Generate a square grid of antennas.

    size (int): Number of antennas in one side of the square grid
    sep (float): seperation in antennas (meters)

    Output: dict of antpos    

    """

    antpos = {}
    for xant in range(side):
        for yant in range(side):
            antpos[xant+yant*side] = np.array([sep*xant, sep*yant, 0.0])

    return antpos

antpos = square_array(6)
Nants = len(antpos)
reds = hera_cal.redcal.get_reds(antpos)
Nbls = float(np.sum([len(s) for s in reds]))
 
Nsim = 256
sigma = 0.01
#true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds, shape=(1,Nsim), 
#                                                            gains = {},
#                                                            gain_scatter=0.1)
#data = {k:v+sigma*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}
#
#num_unibls = np.arange(3,len(reds),1)
#subvisbls = {}
#numbls = {}
#
#
#subsol = {}
#for i in num_unibls:
#    print i
#
#    redset = reds[:i]
#    subvisbls[i] = [r[0] for r in redset]
#    numbls[i] = np.sum([len(s) for s in redset])
#
#    rc = hera_cal.redcal.RedundantCalibrator(redset)
#
#    sd = {k:v for k,v in true_gains.items()}
#    sd.update({k:v for k,v in true_vis.items() if k in subvisbls[i]}) 
#
#    sd = rc.omnical(data, sd)
#    subsol[i] = rc.remove_degen(sd[1], degen_sol=true_gains)
#
#with open('gain_scaling_bls_sqarray_subredcal.cp','w') as fp:
#    cp.dump([antpos, reds, num_unibls, subsol, numbls], fp, protocol=2)

with open('gain_scaling_bls_sqarray_subredcal.cp','r') as fp:
    [antpos, reds, num_unibls, subsol, numbls] = cp.load(fp)

# Plot gain variance vs fraction of baselines

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16, right=0.96, top=0.96, bottom=0.15)

for i in num_unibls:

    g = np.zeros([Nants, Nsim], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant]= subsol[i][(ant,'Jxx')]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    avg = np.mean(g_var)
    ax.errorbar(numbls[i]/Nbls, y= avg, yerr= np.std(g_var),
                fmt='o', color= 'k')

num_unibls = np.asarray(num_unibls, dtype=np.float)
num_redbls = np.asarray([v for v in numbls.values()], dtype=np.float)
ax.plot(num_redbls/Nbls, 0.5*(Nants/num_redbls),'-', 
        color=colors[0])#, label=r'$1/N_{\mathrm{b}}$')
ax.plot(num_redbls/Nbls, 2.5/num_unibls**1.5, '--', 
        color=colors[2])#, label=r'$1/N^{1.5}_{\mathrm{rs}}$')

ax.loglog()
ax.set_xlabel(r'Fraction of all baselines')
ax.set_ylabel(r'$\sigma_{\left| g\right|}^2 $')
ax.grid(which='both', ls='dotted', alpha=0.7)
ax.set_xlim([0.09, 1.3])
ax.set_xticks([0.1,0.5,1])
ax.set_xticklabels([0.1,0.5,1])

plt.show()
