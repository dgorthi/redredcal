import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp
#import hera_sim
import hera_cal
import seaborn as sns

#plt.style.use(['presentation'])
plt.style.use(['paper'])
colors = sns.color_palette('colorblind')

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

NSIM = 256
sigma = 0.3

with open('scale_noise_wgains_hexNum4-26_sig03.cp','r') as fp:
    gains, Nants, nbls_sub = cp.load(fp)

#with open('scale_noise_wchisq_hexNum4-8.cp','r') as fp:
#    gains, mvis, chisq, tg, tv, Nants = cp.load(fp)

# Plot variance of gains
fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16,right=0.95,bottom=0.16,top=0.95)
for hexNum,sol in gains.items():

    antpos = hex_array(hexNum, split_core=False, outriggers=0)

    g = np.zeros([len(antpos), NSIM], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant] = sol['subredcal'][(ant,'Jxx')][0]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    h1 = ax.errorbar(len(antpos), y=np.mean(g_var), yerr=np.std(g_var),
                     fmt='s', ms=5, capsize=5,
                     color=colors[1], label='subredcal')

    g = np.zeros([len(antpos), NSIM], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant] = sol['lowcadcal'][(ant,'Jxx')][0]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    h2 = ax.errorbar(len(antpos), y=np.mean(g_var), yerr=np.std(g_var),
                     fmt='^', ms=5, capsize=5,
                     color=colors[2], label='lowcadcal')

    g = np.zeros([len(antpos), NSIM], dtype=np.complex) 
    for ant in antpos.keys():
        g[ant] = sol['redcal'][(ant,'Jxx')][0]
    g_var = np.abs(np.diag(np.cov(g))) / sigma**2
    h3 = ax.errorbar(len(antpos), y=np.mean(g_var), yerr=np.std(g_var),
                     fmt='o', ms=5, capsize=5,
                     color='k', label='redcal')

Nants = np.asarray(Nants)
ax.loglog(Nants, (Nants**-1.0),  '--', color='k')
ax.loglog(Nants, 0.56/np.log2(Nants), '--', color=colors[2])
ax.loglog(Nants, 0.78/np.log2(Nants), '--', color=colors[1])

ax.legend(handles=[h1,h2,h3], fancybox=True)

ax.grid(ls='dotted', which='both', alpha=0.5)

ax.set_xlabel('Number of antenna elements')
ax.set_ylabel(r'$\sigma_{g}^2 $')
ax.set_ylim([3e-4, 3e-1])

fig.savefig('scaling_antennas.pdf')

plt.show()
