import numpy as np
from hera_cal import redcal
import hera_cal
import matplotlib
import cPickle as cp
import linsolve
import seaborn as sns
import hera_sim

antpos = hera_sim.antpos.hex_array(4, split_core=False, outriggers=0)
reds = redcal.get_reds(antpos)
redbls = redcal.get_pos_reds(antpos)

# Increase SNR by decreasing noise in the measured visibilities
Nsim = 256
Nsim2 = int(np.log(Nsim)/np.log(2))
steps = np.logspace(0, Nsim2, base=2, num=Nsim2+1, dtype=np.int)

true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds, shape=(1,Nsim),
                                                            gains = {},
                                                            gain_scatter=0.01)

sigma = 0.1
data_sub = {k:v+sigma*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}

inttime = np.logspace(-1, 7, num=9, base=2)
sigma_list = 2.5/np.sqrt(inttime) 

lowsol = {}
model_vis = {}
chisq     = {}
rc = redcal.RedundantCalibrator(reds)
for i,sgl in enumerate(sigma_list):
    print "\nSigma Visib: {0:1.3f}, Inttime: {1:2.3f}".format(sgl, inttime[i])
    data_low = {k:v+sgl*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}
    sd = {k:v for k,v in true_gains.items()}
    sd.update({k:v for k,v in true_vis.items()})
    sd = rc.omnical(data_low, sd)
    lowsol[sgl] = rc.remove_degen(sd[1], degen_sol= true_gains)


with open('gain_scaling_snr_lowcadcal.cp','w') as fp:
    cp.dump(lowsol, fp, protocol=2)

