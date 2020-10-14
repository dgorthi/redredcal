import numpy as np
import cPickle as cp
import hera_sim
from hera_cal import redcal

Nsim = 256
sigma_list = np.logspace(-1, 1, num=11)
lowsol = {}

hexNums = np.arange(4,10)
for hexNum in hexNums:
    antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
    reds = redcal.get_reds(antpos)

    print hexNum, len(antpos)
    
    true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds,
                                                                shape=(1,Nsim),
                                                                gains={},
                                                                gain_scatter=0.01)
    rc = redcal.RedundantCalibrator(reds)
    lowsol[len(antpos)] = {}

    for sigma in sigma_list:
        print sigma
        data = {k:v+sigma*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}
        sd = {k:v for k,v in true_gains.items()}
        sd.update({k:v for k,v in true_vis.items()})
        sd = rc.lincal(data, sd)
        lowsol[len(antpos)][sigma] = rc.remove_degen(sd[1], degen_sol=true_gains)

    with open('gain_scaling_snr_ants_lowcadcal_lincal.cp','w') as fp:
        cp.dump(lowsol, fp, protocol=2)  
