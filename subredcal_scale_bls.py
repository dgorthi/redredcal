import numpy as np
import hera_cal
import matplotlib
import cPickle as cp
import linsolve
import seaborn as sns
#import hera_sim
from copy import deepcopy

def hex_array(hex_num, sep=14.6, split_core=True, outriggers=2):
    """
    Build a hexagonal array configuration, nominally matching HERA's ideal configuration.

    Args:
        hex_num (int): the hexagon (radial) number of the core configuration.
            Number of core antennas returned is 3N^2 - 3N + 1.
        sep (float): the separation between hexagonal grid points (meters).
        split_core (bool): fractures the hexagonal core into tridrents that subdivide
            a hexagonal grid. Loses N antennas, so the number of core antennas returned
            is 3N^2 - 4N + 1.
        outriggers (int): adds R extra rings of outriggers around the core that tile
            with the core to produce a fully-sampled UV plane.  The first ring
            corresponds to the exterior of a hexNum=3 hexagon. Adds 3R^2 + 9R antennas.
    Returns:
        dict: a dictionary of antenna numbers and positions.
            Positions are x,y,z in topocentric coordinates, in meters.
    """
    # Main Hex
    positions = []
    for row in range(
            hex_num - 1, -hex_num + split_core, -1
    ):  # the + split_core deletes a row
        for col in range(0, 2 * hex_num - abs(row) - 1):
            x_pos = sep * ((-(2 * hex_num - abs(row)) + 2) / 2.0 + col)
            y_pos = row * sep * 3 ** 0.5 / 2
            positions.append([x_pos, y_pos, 0])

    # unit vectors
    up_right = sep * np.asarray([0.5, 3 ** 0.5 / 2, 0])
    up_left = sep * np.asarray([-0.5, 3 ** 0.5 / 2, 0])

    # Split the core into 3 pieces
    if split_core:
        new_pos = []
        for i, pos in enumerate(positions):
            theta = np.arctan2(pos[1], pos[0])
            if pos[0] == 0 and pos[1] == 0:
                new_pos.append(pos)
            elif -np.pi / 3 < theta < np.pi / 3:
                new_pos.append(np.asarray(pos) + (up_right + up_left) / 3)
            elif np.pi / 3 <= theta < np.pi:
                new_pos.append(np.asarray(pos) + up_left - (up_right + up_left) / 3)
            else:
                new_pos.append(pos)
        positions = new_pos

    # Add outriggers
    if outriggers:
        exterior_hex_num = outriggers + 2
        for row in range(exterior_hex_num - 1, -exterior_hex_num, -1):
            for col in range(2 * exterior_hex_num - abs(row) - 1):
                x_pos = (
                        ((-(2 * exterior_hex_num - abs(row)) + 2) / 2.0 + col)
                        * sep
                        * (hex_num - 1)
                )
                y_pos = row * sep * (hex_num - 1) * 3 ** 0.5 / 2
                theta = np.arctan2(y_pos, x_pos)
                if (x_pos ** 2 + y_pos ** 2) ** 0.5 > sep * (hex_num + 1):
                    # These specific displacements of the outrigger sectors are designed specifically
                    # for redundant calibratability and "complete" uv-coverage, but also to avoid
                    # specific obstacles on the HERA site (e.g. a road to a MeerKAT antenna).
                    if 0 < theta <= 2 * np.pi / 3 + 0.01:
                        positions.append(
                            np.asarray([x_pos, y_pos, 0]) - 4 * (up_right + up_left) / 3
                        )
                    elif 0 >= theta > -2 * np.pi / 3:
                        positions.append(
                            np.asarray([x_pos, y_pos, 0]) - 2 * (up_right + up_left) / 3
                        )
                    else:
                        positions.append(
                            np.asarray([x_pos, y_pos, 0]) - 3 * (up_right + up_left) / 3
                        )

    return {i: pos for i, pos in enumerate(np.array(positions))}

# Load array
antpos = hera_sim.antpos.hex_array(4, split_core=False, outriggers=0)
Nants = len(antpos)
reds = hera_cal.redcal.get_reds(antpos)

Nsim = 256
sigma = 0.001
true_gains, true_vis, true_data = hera_sim.vis.sim_red_data(reds, shape=(1,Nsim), 
                                                            gains = {},
                                                            gain_scatter=0.0)
data = {k:deepcopy(v)+sigma*hera_sim.noise.white_noise((1,Nsim)) for k,v in true_data.items()}


# -------------------------------------------------
# Increment baselines in order of baseline length
# -------------------------------------------------

print 'Increment baselines in order'

uni_bl_types = np.arange(len(reds))
num_uni_bls  = np.arange(1,len(reds)+1,1)

ants_per_unibl = {} # The number of antennas a baseline type touches
bls_per_ant  = {} # Number of baselines per antenna
subvisbls = {}
numbls = {}

for num,bltype in zip(num_uni_bls, uni_bl_types):
    ants_per_unibl[bltype] = len(np.unique([[a[0],a[1]] for a in reds[bltype]]))
    
    if num<3: continue
    redset = [reds[x] for x in uni_bl_types[:num]]
    subvisbls[bltype] = [r[0] for r in redset]
    numbls[num]    = np.sum([len(s) for s in redset])

    x = []
    for ant in range(Nants):
        x.append(len([bl for subbl in redset for bl in subbl if ant in bl]))
    bls_per_ant[num] = np.mean(x)

    print '{:2d}\t{:2d}\t{:.3f}\t{:3d}'.format(num, numbls[num], bls_per_ant[num],
                                              ants_per_unibl[bltype])

subsol = {}
for num, bltype in zip(num_uni_bls, uni_bl_types):
    if num<3: continue
    print '{:2d}\t{:2d}'.format(num, bltype)
    redset = [reds[x] for x in uni_bl_types[:num]]
    rc = hera_cal.redcal.RedundantCalibrator(redset)

    sd = deepcopy({k:v for k,v in true_vis.items() if k in subvisbls[bltype]})
    sd.update(deepcopy(true_gains))

    sd = rc.lincal(data, sd)
    subsol[bltype] = rc.remove_degen(sd[1], degen_sol = true_gains)

with open('gain_scaling_bls_subredcal.cp', 'w') as fp:
    cp.dump([subsol, uni_bl_types, subvisbls, numbls, ants_per_unibl, bls_per_ant], 
            fp, protocol=2)

## -------------------------------------------------
## Increment baselines in random order
## -------------------------------------------------
#
#print 'Increment baselines randomly'
#
#uni_bl_types = np.arange(3, len(reds), 1)
#np.random.shuffle(uni_bl_types)
#uni_bl_types = np.append(np.arange(3), uni_bl_types)
#
#ants_per_unibl = {} # The number of antennas a baseline type touches
#subvisbls = {}
#numbls = {}
#
#for num,bltype in zip(num_uni_bls, uni_bl_types):
#    ants_per_unibl[bltype] = len(np.unique([[a[0],a[1]] for a in reds[bltype]]))
#
#    if num<3: continue
#    redset = [reds[x] for x in uni_bl_types[:num]]
#    subvisbls[bltype] = [r[0] for r in redset]
#    numbls[num]    = np.sum([len(s) for s in redset])
#
#    print '{:2d}\t{:2d}\t{:3d}\t{:3d}'.format(num, bltype, numbls[num], 
#                                              ants_per_unibl[bltype])
#
#subsol = {}
#for num, bltype in zip(num_uni_bls, uni_bl_types):
#    if num<3: continue
#    print '{:2d}\t{:2d}'.format(num, bltype)
#    redset = [reds[x] for x in uni_bl_types[:num]]
#    rc = hera_cal.redcal.RedundantCalibrator(redset)
#
#    sd = deepcopy({k:v for k,v in true_vis.items() if k in subvisbls[bltype]})
#    sd.update(deepcopy(true_gains))
#
#    sd = rc.omnical(data, sd)
#    subsol[bltype] = rc.remove_degen(sd[1], degen_sol = true_gains)
#
#with open('gain_scaling_bls_random_subredcal.cp', 'w') as fp:
#    cp.dump([subsol, uni_bl_types, subvisbls, numbls, ants_per_unibl], fp, protocol=2)
#


#num_redbls = np.arange(3,63,1)
#np.random.shuffle(num_redbls) #inplace operation
#num_redbls = np.append(np.asarray([0,1,2]), num_redbls) #for degeneracy criterion
#
#ants_per_bl = {}
#subvisbls = {}
#numbls = {}
#
#for i in num_redbls:
##for i in range(3,len(num_redbls)):
#    redset = reds[:i]
#    #redset = [reds[x] for x in num_redbls[:i]]
#    subvisbls[i] = [r[0] for r in redset]
#    numbls[i] = len([pair for subbl in redset for pair in subbl])
#
#    #ants_per_bl[i] = len(np.unique(np.append([a[0] for a in reds[num_redbls[i]]],
#    #                                           [a[1] for a in reds[num_redbls[i]]])))
#
#    #print i, num_redbls[i], numbls[i], ants_per_bl[i] 
#
##ants_per_bl_arr = np.append([Nants, Nants, Nants],[v for v in ants_per_bl.values()])
#
#subsol = {}
#for i in num_redbls:
##for i in range(3,len(num_redbls)):
#    print i
#    redset = reds[:i]
#    #redset = [reds[x] for x in num_redbls[:i]]
#    rc = hera_cal.redcal.RedundantCalibrator(redset)
#    sd = {k:v for k,v in true_gains.items()}
#    sd.update({k:v for k,v in true_vis.items() if k in subvisbls[i]}) 
#    sd = rc.omnical(data, sd)
#    subsol[numbls[i]] = rc.remove_degen(sd[1], degen_sol=true_gains)
#    #subsol[i] = rc.remove_degen(sd[1], degen_sol= true_gains)
#
#with open('gain_scaling_bls_legacy_subredcal.cp','w') as fp:
#    cp.dump([subsol, num_redbls], fp, protocol=2)

