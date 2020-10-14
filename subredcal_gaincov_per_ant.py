import numpy as np
import matplotlib.pyplot as plt
from hera_cal import redcal
import hera_cal
from matplotlib import colors

#plt.style.use(['dark_background','presentation'])
plt.style.use(['paper'])
## Cycler colors because I dont know how to extract them from the style file
#colors = ['#0173b2','#de8f05','#029e73','#d55e00','#cc78bc',
#          '#ca9161','#fbafe4','#949494','#ece133','#56b4e9']
line_style=[':','--','-']
marker = ['o','^','s']

def genA(redbls):
    # number of parameters = num_ants + num_unique baselines
    N = Nants + len(redbls)
    
    # number of measurements = total number of baselines
    allbls = [bl for red in redbls for bl in red]
    M = len(allbls)
    
    A = np.zeros([M,N],dtype=np.complex)
    i = 0
    for bl,reds in enumerate(redbls):
        for pair in reds:
            A[i,pair[0]] = (1+1j)
            A[i,pair[1]] = (1-1j)
            A[i,Nants+bl] = 1+1j
            i += 1
    return np.matrix(A)

# Setup array
ants = np.loadtxt('antenna_positions_37.dat')
antpos = {k:v for k,v in zip(range(len(ants)), ants)}
Nants = len(antpos)

redbls = redcal.get_pos_reds(antpos)
reds = redcal.get_reds(antpos)

# Redcal
A  = genA(redbls)
Mr = np.dot(np.real(A.T), np.real(A))
Mi = np.dot(np.imag(A.T), np.imag(A))
covr = np.linalg.pinv(Mr)[:Nants,:Nants]
covi = np.linalg.pinv(Mi)[:Nants,:Nants]

# Subredcal
At = genA(redbls[0:3])
Mrt = np.dot(np.real(At.T), np.real(At))
Mit = np.dot(np.imag(At.T), np.imag(At))
covtr = np.linalg.pinv(Mrt)[:Nants,:Nants]
covti = np.linalg.pinv(Mit)[:Nants,:Nants]

# Subredcal v2
Att = genA(redbls[0:17])
Mrtt = np.dot(np.real(Att.T), np.real(Att))
Mitt = np.dot(np.imag(Att.T), np.imag(Att))
covttr = np.linalg.pinv(Mrtt)[:Nants,:Nants]
covtti = np.linalg.pinv(Mitt)[:Nants,:Nants]

fig,ax = plt.subplots(3,2,figsize=(6.4,7.2))
fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05,
                    hspace=0.1, wspace=0)
cmap = 'magma_r'

# Real
im = ax[0][0].imshow(covr, cmap=cmap)
cb = fig.colorbar(im, ax=ax[0][0])
cb.ax.tick_params(labelsize=12)

im = ax[1][0].imshow(covtr,  cmap=cmap)
cb = fig.colorbar(im, ax=ax[1][0])
cb.ax.tick_params(labelsize=12)

im = ax[2][0].imshow(covttr, cmap=cmap)
cb = fig.colorbar(im, ax=ax[2][0])
cb.ax.tick_params(labelsize=12)

# Imag
im = ax[0][1].imshow(covi, cmap=cmap)
cb = fig.colorbar(im, ax=ax[0][1])
cb.ax.tick_params(labelsize=12)

im = ax[1][1].imshow(covti, cmap=cmap)
cb = fig.colorbar(im, ax=ax[1][1])
cb.ax.tick_params(labelsize=12)

im = ax[2][1].imshow(covtti, cmap=cmap)
cb = fig.colorbar(im, ax=ax[2][1])
cb.ax.tick_params(labelsize=12)

for row in ax:
    for col in row:
        col.set_aspect('equal')
        col.tick_params(left=False, right=False, 
                        top=False, bottom=False,
                        labelleft=False, labelright=False,
                        labeltop=False, labelbottom=False)

# Labels
ax[0][0].set_title('Log-Amplitude')
ax[0][1].set_title('Phase')
ax[0][0].set_ylabel('(a)', rotation=0, labelpad=20)
ax[1][0].set_ylabel('(b)', rotation=0, labelpad=20)
ax[2][0].set_ylabel('(c)', rotation=0, labelpad=20)

fig.savefig('covariance.pdf',bbox_inches='tight')

# ----------------------------------
# Gain covariance for one antenna
# ----------------------------------
x = [v[0] for v in antpos.values()]
y = [v[1] for v in antpos.values()]

fig,ax = plt.subplots(2,2)
fig.subplots_adjust(left=0.02,right=0.86,top=0.98,bottom=0.03,hspace=0.08,wspace=0.04)
cmap = 'magma_r'
maxc = np.max(np.abs(covtr+1j*covti))
minc = np.min(np.abs(covtr+1j*covti))
#norm = colors.Normalize(vmin=minc, vmax=maxc)
#norm = colors.SymLogNorm(linthresh=0.01, linscale=1, vmin=minc, vmax=maxc)
norm = colors.LogNorm(vmin=0.025, vmax=maxc)

for i,ant in enumerate([0,1,5,18]):
    print (i//2, i%2)
    c = np.abs(np.asarray(covtr[ant])[0] + 1j*np.asarray(covti[ant])[0])
    # see https://en.wikipedia.org/wiki/Complex_random_vector#Cross-covariance_matrix_and_pseudo-cross-covariance_matrix
    print (c)
    im = ax[i//2][i%2].scatter(x,y,c=c, marker='o', edgecolors='k', 
                               norm = norm,
                               #norm=matplotlib.colors.LogNorm(vmin=0.04, vmax=0.6), 
                               cmap=cmap, s=550)
    ax[i//2][i%2].plot(x[ant],y[ant],'*w',markersize=8)

for row in ax:
    for col in row:
        col.set_aspect('equal')
        col.set_xlim([-110.2,-6.6])
        col.set_ylim([-121.8,-29.9])
        col.tick_params(left=False, right=False, 
                        top=False, bottom=False,
                        labelleft=False, labelright=False,
                        labeltop=False, labelbottom=False)

cax = fig.add_axes([0.87,0.03,0.015,0.95]) 
cbar = fig.colorbar(im, cax=cax, ticks=[0.03, 0.05, 0.1, 0.25, 0.5], format="%.2f")
cbar.ax.set_yticklabels(['%.2f'%s for s in [0.03, 0.05, 0.1, 0.25, 0.5]])
fig.savefig('gaincov.pdf')

plt.show()
