import numpy as np
from hera_cal import redcal
import linsolve
import cPickle as cp
from pyuvdata.utils import jnum2str, jstr2num
import logging
import hera_sim
import hera_cal
from copy import deepcopy

logging.basicConfig(level= logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')

def optimize_red_vis(data, gains, subreds):
    eqs = {}; const = {}
    for (i,j,pol) in subreds:
        eqs['g_%d * g_%d_ * V'%(i,j)] = data[(i,j,pol)]

    for i in range(Nants):
        const['g_%d'%(i)] = gains[(i,'Jxx')]    

    lps = linsolve.LinearSolver(eqs, **const)
    X = lps.solve()
    return X['V']

Nsim = 256
hexNumrange = np.arange(4,8,1)
model_vis = {}
gains = {}
Nants_list = []
chi = {}
tg = {}
tv = {}

## Tmax sets the time over which the baselines are computed for lowcadcal
## or the time for which the gains are averaged under subredcal. Here I have
## translated Tmax into a number of simulations unit.
## Each simulation ~ 10sec
## 40   sec === 4 sim
## 5.3  min === 32 sim
## 10.6 min === 64 sim
Tsim = 1

logger = logging.getLogger('scaling_ants_sim')
logger.setLevel(logging.DEBUG)

for hexNum in hexNumrange:

    logging.info('hexNum: %d, Generating hex array..'%hexNum)
    antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
    logging.info('hexNum: %d, Computing redundant baselines..'%hexNum)
    reds = hera_cal.redcal.get_reds(antpos)
    
    Nants = len(antpos)
    Nants_list.append(Nants)
    Nvars = Nants + len(reds)
    Nbls  = np.sum([len(s) for s in reds])

    model_vis[hexNum] = {}
    chi[hexNum] = {}
    gains[hexNum] = {}
    tg[hexNum] = {}
    tv[hexNum] = {}

    logging.info('hexNum: %d,  Nants: %d,  Nbls: %d' %(hexNum, Nants, Nbls))

    # ---------------------------------    
    #  Simulate gains and visibilities
    # ---------------------------------
    tg[hexNum], tv[hexNum], td = hera_sim.vis.sim_red_data(reds, shape=(1,Nsim),  
                                                           gain_scatter=0.1)

    #tg[hexNum] = {k: np.repeat(v, Nsim, axis=1) for k,v in tg[hexNum].items()}   
    #tv[hexNum] = {k: np.repeat(v, Nsim, axis=1) for k,v in tv[hexNum].items()}
    #td = {k: np.repeat(v, Nsim, axis=1) for k,v in td.items()}
     
    sigma = 0.3
    data = {k:v+ sigma*hera_sim.noise.white_noise((1, Nsim)) for k,v in td.items()}
    
    # ----------------------------------
    #             Redcal
    # ----------------------------------
    logging.info('redcal: solving..')
    rc = redcal.RedundantCalibrator(reds)
    sd = deepcopy({k:v for k,v in tv[hexNum].items()})
    sd.update(deepcopy({k:v for k,v in tg[hexNum].items()}))

    logging.info('redcal: Omnical..')
    sd = rc.omnical(data, sd, maxiter=1024)

    logging.info('redcal: Remove degeneracy')
    sol = rc.remove_degen(sd[1], degen_sol=tg[hexNum])
    gains[hexNum]['redcal'] = deepcopy(sol)
     
    logging.info('redcal: Compute chisq')
    chisq = 0
    for k in tv[hexNum].keys():
        subbls = [bls for bls in reds if k in bls][0]
        for (i, j, pol) in subbls:
            chisq += np.abs((sol[i,'Jxx']*np.conj(sol[j,'Jxx'])*sol[k]) - data[(i,j,pol)])**2
    
    chi[hexNum]['redcal'] = chisq/((Nbls-Nvars-2)*sigma**2)
    
    # -----------------------------------
    #            Subredcal
    # -----------------------------------

    # Compute the number of baselines that make the 
    # correlator NlogN

    logging.info('subredcal: Starting!')
    NlogN = int(3* Nants * np.log2(2*Nants))
    N = 3
    subreds = reds[0:N]
    Ncorr = np.sum([len(s) for s in subreds])
    while (Ncorr < NlogN):
        N = N+1
        subreds = reds[0:N]
        Ncorr = np.sum([len(s) for s in subreds])
    
    N = N-1
    subreds = reds[0:N]
    subvisbls = [r[0] for r in subreds]
    Ncorr = np.sum([len(s) for s in subreds])
    logging.info('subredcal: Unique baseline types included = %d' % N)
    logging.info('subredcal: Number of baselines = %d' % Ncorr)
    logging.info('subredcal: assert(baselines<NlogN): %s' % (Ncorr < NlogN))
    
    rc = redcal.RedundantCalibrator(subreds)

    sd = deepcopy({k:v for k,v in tg[hexNum].items()})
    sd.update(deepcopy({k:v for k,v in tv[hexNum].items() if k in subvisbls}))

    logging.info('subredcal: solving..')
    sd = rc.omnical(data, sd)

    subsol = rc.remove_degen(sd[1], degen_sol=tg[hexNum])
    gains[hexNum]['subredcal'] = deepcopy(subsol)
    
    logging.info('subredcal: Compute model visibilities')

    # Compute all visibilities using averaged gains
    model_vis_sub = {}
    for k in tv[hexNum].keys():
        if k in subsol.keys():
            model_vis_sub[k] = subsol[k]
        else:
            subbls = [bls for bls in reds if k in bls][0]
            model_vis_sub[k] = optimize_red_vis(data, subsol, subbls)

    model_vis[hexNum]['subredcal'] = deepcopy(model_vis_sub)
    
    logging.info('subredcal: Compute chisq')
    chisq = 0
    for k in tv[hexNum].keys():
        subbls = [bls for bls in reds if k in bls][0]
        for (i,j,pol) in subbls:
            chisq += np.abs((subsol[i,'Jxx']*np.conj(subsol[j,'Jxx'])*model_vis_sub[k]) - data[(i,j,pol)])**2
    
    chi[hexNum]['subredcal'] = chisq/((Nbls-Nvars-2)*sigma**2)
    
    chisq = 0
    for k in subvisbls:
        subbls = [bls for bls in reds if k in bls][0]
        for (i,j,pol) in subbls:
            chisq += np.abs((subsol[i,'Jxx']*np.conj(subsol[j,'Jxx'])*subsol[k]) - data[(i,j,pol)])**2
    
    chi[hexNum]['subredcal_sub'] = chisq/((Ncorr-Nants-N-2)*sigma**2)
    

    # --------------------------------------
    #               Lowcadcal
    # --------------------------------------
    
    # Generate noisier data to reflect the lower integration time
    sigma_low = sigma * np.sqrt(1.0*Nbls/Ncorr)
    logging.info('lowcadcal: sigma=%1.4f'%sigma_low)
    data_low = {k:v+sigma_low*hera_sim.noise.white_noise((1,Nsim)) for k,v in td.items()}
    
    logging.info('lowcadcal: solving..')
    rc = redcal.RedundantCalibrator(reds)
 
    sd = deepcopy({k:v for k,v in tg[hexNum].items()})
    sd.update(deepcopy({k:v for k,v in tv[hexNum].items()}))

    logging.info('lowcadcal: Omnical..')
    sd = rc.omnical(data_low, sd, maxiter=1024)

    logging.info('lowcadcal: Remove Degen')
    lowsol = rc.remove_degen(sd[1], degen_sol= tg[hexNum])

    gains[hexNum]['lowcadcal'] = deepcopy(lowsol)
    
    # Compute model visibilities
    logging.info('lowcadcal: Compute model vis')
    model_vis_low = {} 
    for k in tv[hexNum].keys():
        subbls = [bls for bls in reds if k in bls][0]
        model_vis_low[k] = optimize_red_vis(data, lowsol, subbls)

    model_vis[hexNum]['lowcadcal'] = deepcopy(model_vis_low)
    
    logging.info('lowcadcal: Compute chisq')
    chisq = 0
    for k in tv[hexNum].keys():
        subbls = [bls for bls in reds if k in bls][0]
        for (i,j,pol) in subbls:
            chisq += np.abs((lowsol[i,'Jxx']*np.conj(lowsol[j,'Jxx'])*model_vis_low[k]) - data[(i,j,pol)])**2
    
    chi[hexNum]['lowcadcal'] = chisq/((Nbls-Nvars-2)*sigma**2)

    chisq = 0
    for k in tv[hexNum].keys():
        subbls = [bls for bls in reds if k in bls][0]
        for (i,j,pol) in subbls:
            chisq += np.abs((lowsol[i,'Jxx']*np.conj(lowsol[j,'Jxx'])*lowsol[k]) - data_low[(i,j,pol)])**2
    
    chi[hexNum]['lowcadcal_low'] = chisq/((Nbls-Nvars-2)*sigma_low**2)


    with open('scale_noise_wchisq_hexNum4-8.cp','w') as fp:
        cp.dump([gains, model_vis, chi, tg, tv, Nants_list], fp, protocol=2)
