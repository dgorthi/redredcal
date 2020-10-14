import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from matplotlib.patches import Polygon

plt.style.use(['paper'])
colors = sns.color_palette('colorblind')

c_sub = colors[1]
c_low = colors[2]
c_red = 'k'
c_box = 'silver'

#with open('compare_gainvar_chisq_hexNums4-40_final.pkl','rb') as fp:
#    gains, chi, Nubls = pickle.load(fp)

with open('compare_gainvar_boxplot_hexNums4-60.pkl', 'rb') as fp:
    gvar = pickle.load(fp)

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.16,right=0.95,bottom=0.16,top=0.95)

Nants_list = np.sort(list(gvar.keys()))

Nants_plot_list = [37, 61, 91, 127, 169, 217, 271, 331, 397, 469, 547, 
                   631, 721, 817, 919, 1027, 1261, 1519, 
                   1801, 2107, 2611, 3169, 3571, 3997, 4681, 7351, 10621]

w = 0.03
width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)

for nants in Nants_plot_list:
    
    h3 = ax.boxplot(gvar[nants]['redcal'], positions=[nants], widths=width(nants,w),
                    showfliers=False, patch_artist=True,
                    showmeans=True,
                    meanprops=dict(marker='o', color=c_red, markersize=3,
                                   markeredgecolor=c_red, markerfacecolor=c_red),
                    medianprops=dict(linestyle='-', color=c_box))
    plt.setp(h3['boxes'],    color=c_box) 
    plt.setp(h3['whiskers'], color=c_box) 
    plt.setp(h3['caps'],     color=c_box)

    box = h3['boxes'][0]
    box.set_facecolor(c_box)

    caps = h3['caps']
    for c in caps:
        c.set(xdata=c.get_xdata()+ (-width(nants,0.01),width(nants,0.01)))

    h1 = ax.boxplot(gvar[nants]['subredcal'], positions=[nants], 
                    widths=width(nants,w), showfliers=False,
                    showmeans=True,
                    meanprops=dict(marker='s', color=c_sub, markersize=3,
                                   markeredgecolor=c_sub, markerfacecolor=c_sub),
                    medianprops=dict(linestyle='-', color=c_box))
    plt.setp(h1['boxes'],    color=c_box) 
    plt.setp(h1['whiskers'], color=c_box) 
    plt.setp(h1['caps'],     color=c_box) 

    box = h1['boxes'][0]
    box_coords = np.transpose(np.stack([box.get_xdata(), box.get_ydata()]))
    ax.add_patch(Polygon(box_coords, facecolor=c_box))

    caps = h1['caps']
    for c in caps:
        c.set(xdata=c.get_xdata()+ (-width(nants,0.01),width(nants,0.01)))


    h2 = ax.boxplot(gvar[nants]['lowcadcal'], positions=[nants+0.5*width(nants,w)], 
                    widths=width(nants,w), showfliers=False,
                    showmeans=True,
                    meanprops=dict(marker='^', color=c_low, markersize=3,
                                   markeredgecolor=c_low, markerfacecolor=c_low),
                    medianprops=dict(linestyle='-', color=c_box))
    plt.setp(h2['boxes'],    color=c_box)
    plt.setp(h2['whiskers'], color=c_box)
    plt.setp(h2['caps'],     color=c_box)

    box = h2['boxes'][0]
    box_coords = np.transpose(np.stack([box.get_xdata(), box.get_ydata()]))
    ax.add_patch(Polygon(box_coords, facecolor=c_box))
    
    caps = h2['caps']
    for c in caps:
        c.set(xdata=c.get_xdata()+ (-width(nants,0.01),width(nants,0.01)))


    #h1 = ax.errorbar(Nants, y= gains[Nants]['subredcal']['avg'], 
    #                 yerr= 0.5*gains[Nants]['subredcal']['std'],
    #                 fmt='s', ms=5, capsize=5,
    #                 color=colors[1], label='subredcal')

    #h2 = ax.errorbar(Nants, y= gains[Nants]['lowcadcal']['avg'], 
    #                 yerr= 0.5*gains[Nants]['lowcadcal']['std'],
    #                 fmt='^', ms=5, capsize=5,
    #                 color=colors[2], label='lowcadcal')

    #h3 = ax.errorbar(Nants, y= gains[Nants]['redcal']['avg'], 
    #                 yerr= 0.5*gains[Nants]['redcal']['std'],
    #                 fmt='o', ms=5, capsize=5,
    #                 color='k', label='redcal')

ax.loglog(Nants_list, (Nants_list**-1.0),  '--', color=c_red, lw=1)
ax.loglog(Nants_list, 0.50/np.log2(Nants_list), '--', color=c_low, lw=1)
ax.loglog(Nants_list, 0.75/np.log2(Nants_list), '--', color=c_sub, lw=1)

#ax.legend(handles=[h1,h2,h3], fancybox=True)
ax.plot(0,0,'s', markersize=5, color=c_sub, label='subredcal')
ax.plot(0,0,'^', markersize=5, color=c_low, label='lowcadcal')
ax.plot(0,0,'o', markersize=5, color=c_red, label='redcal')
ax.legend(fancybox=True)

ax.grid(ls='dotted', which='both', alpha=0.5)

ax.set_xlabel('Number of antenna elements')
ax.set_ylabel('$\sigma_{g}^2 $')
ax.set_ylim([2e-5, 3e-1])

fig.savefig('scale_ants_gainvar.pdf')

plt.show()
