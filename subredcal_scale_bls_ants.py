import numpy as np
from hera_cal import redcal
import cPickle as cp
import hera_sim

Nsim = 256
hexNums = np.asarray([4,5,7,9,11,13])
subsol = {}

def get_subreds(num_bls, reds):
    """ Get the list of redundant sets that have 
        number of baseslines <= num_bls
    """
    N=0; i=0;
    for num, subbls in enumerate(reds):
        if N >= num_bls:
           return num
        N += len(subbls)

fraction_bls = [0.12, 0.17, 0.2, 0.24, 0.27, 
                0.31, 0.35, 0.4, 0.45, 0.53,
                0.62, 0.69, 0.8, 0.87, 0.99]

for hexNum in hexNums:
    antpos = hera_sim.antpos.hex_array(hexNum, split_core=False, outriggers=0)
    reds = redcal.get_reds(antpos)
    Nbls = len([pair for subbls in reds for pair in subbls])

    print hexNum, len(antpos), Nbls

    true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds,
                                                                shape=(1,Nsim),
                                                                gains={},
                                                                gain_scatter=0.1)
    sigma = 0.3
    data = {k:v+ sigma*hera_sim.noise.white_noise((1, Nsim)) for k,v in true_data.items()}

    subsol[hexNum] = {}

    for f in fraction_bls:
        N = get_subreds(f*Nbls, reds)
        print "{:0.2f}\t{:3d}".format(f, N)
        
        redset = reds[0:N]
        subvisbls = [r[0] for r in redset]
        rc = redcal.RedundantCalibrator(redset)
        sd = {k:v for k,v in true_gains.items()}
        sd.update({k:v for k,v in true_vis.items() if k in subvisbls})

        sd = rc.omnical(data, sd)
        subsol[hexNum][f] = rc.remove_degen(sd[1], degen_sol=true_gains)

    with open('gain_scaling_bls_ants_subredcal_omnical.cp','w') as fp:
        cp.dump(subsol, fp, protocol=2)


    

    
