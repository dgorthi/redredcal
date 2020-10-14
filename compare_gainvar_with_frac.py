import numpy as np
import matplotlib.pyplot as plt
from hera_cal import redcal
import hera_sim
import hera_cal
import cPickle as cp
from copy import deepcopy

def get_subreds(num_bls, reds):
    """ Get the list of redundant sets that have 
        number of baseslines <= num_bls
    """
    N=0; i=0;
    for num, subbls in enumerate(reds):
        if N >= num_bls:
           return num
        N += len(subbls)

def calibrate(data, gains):
    calib_data = {}
    for (i,j,pol) in data.keys():
        calib_data[(i,j,pol)] = data[(i,j,pol)]/(gains[(i,'Jxx')]*np.conj(gains[(j,'Jxx')]))
    return calib_data

def average_red_vis(data, gains, subreds):
    sub_bl_visib = calibrate(data, gains)
    vis = 0
    for subbl in subreds:
        vis += sub_bl_visib[subbl]/len(subreds)
    return vis

hexrange = np.asarray([4,5,7,9])
Nsim = 256
sigma = 0.01
fraction_bls = np.arange(0.1, 0.6, step=0.1)
subsol = {}
lowsol = {}
chisub = {}
chilow = {}

for hexNum in hexrange:
    print 'Processing hexNum {:d}\n'.format(hexNum)

    antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
    Nants = len(antpos)

    reds = redcal.get_reds(antpos)
    Nbls = len([pair for subbls in reds for pair in subbls])

    true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds, shape=(1, Nsim), 
                                                                gain_scatter=0.1)
    data = {k:v+ sigma*hera_sim.noise.white_noise((1, Nsim)) for k,v in true_data.items()}

    sigma_range = []

    # -------------
    #   SUBREDCAL
    # -------------
    
    subsol[hexNum] = {}
    chisub[hexNum] = {}

    for f in fraction_bls:
        N = get_subreds(f*Nbls, reds)
        print "{:0.2f}\t{:3d}\t".format(f, N)

        redset = reds[0:N]
        subvisbls = [r[0] for r in redset]
        Nvars = Nants + len(redset)

        rc = redcal.RedundantCalibrator(redset)
        sd = deepcopy(true_gains)
        sd.update(deepcopy({k:v for k,v in true_vis.items() if k in subvisbls}))

        sd = rc.omnical(data, sd)
        sol = rc.remove_degen(sd[1], degen_sol=true_gains)

        subsol[hexNum][f] = deepcopy(sol)

        # Compute model visibilities
        model_vis_sub = {}
        for k in true_vis.keys():
            if k in sol.keys():
                model_vis_sub[k] = sol[k]
            else:
                subbls = [bls for bls in reds if k in bls][0]
                model_vis_sub[k] = average_red_vis(data, sol, subbls)

        # Compute chisq of fit
        chisq = 0
        for k in true_vis.keys():
            subbls = [bls for bls in reds if k in bls][0]
            for (i,j,pol) in subbls:
                chisq += np.abs((sol[i,'Jxx']*np.conj(sol[j,'Jxx'])*model_vis_sub[k]) - data[(i,j,pol)])**2
        
        chisub[hexNum][f] = chisq/((Nbls-Nvars-2)*sigma**2)

        # For lowcadcal:
        Ncorr = len([pair for subbls in redset for pair in subbls])
        sigma_range.append(sigma * np.sqrt(1.0*Nbls/Ncorr))
        

    # --------------
    #   LOWCADCAL
    # --------------

    lowsol[hexNum] = {}
    chilow[hexNum] = {}
    rc = redcal.RedundantCalibrator(reds)
    
    for sig in sigma_range:
        print "{:1.3f}".format((sigma/sig)**2)
        dlow = {k:v+sig*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}
        
        sd = deepcopy(true_gains)
        sd.update(deepcopy(true_vis))
        
        sd = rc.omnical(dlow, sd)
        sol = rc.remove_degen(sd[1], degen_sol=true_gains)
        lowsol[hexNum]['%.2f'%(sigma/sig)**2] = deepcopy(sol)

        # Compute visibilities
        model_vis_low = {} 
        for k in true_vis.keys():
            subbls = [bls for bls in reds if k in bls][0]
            model_vis_low[k] = average_red_vis(data, sol, subbls)
        
        # Compute Chisq
        chisq = 0
        for k in true_vis.keys():
            subbls = [bls for bls in reds if k in bls][0]
            for (i,j,pol) in subbls:
                chisq += np.abs((sol[i,'Jxx']*np.conj(sol[j,'Jxx'])*model_vis_low[k]) - data[(i,j,pol)])**2
        
        chilow[hexNum]['%.2f'%(sigma/sig)**2] = chisq/((Nbls-Nvars-2)*sigma**2)

    with open('lowcadcal_vs_subredcal_ants.cp','w') as fp:
         cp.dump([subsol, chisub, lowsol, chilow], fp, protocol=2)

