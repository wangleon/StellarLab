#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from stella.parameter.metal import feh_to_z

fig = plt.figure(figsize=(10,4), dpi=150)
ax1 = fig.add_axes([0.07,0.15,0.40,0.80])
ax2 = fig.add_axes([0.50,0.15,0.40,0.80], projection='3d')
ax3 = fig.add_axes([0.90,0.15,0.02,0.80])

fe0, fe1, dfe = -3.0, 1.0, 0.10
af0, af1, daf =  0.0, 0.6, 0.02
feh_lst   = np.arange(fe0, fe1+1e-6, dfe)
alpha_lst = np.arange(af0, af1+1e-6, daf)
feh, alpha = np.meshgrid(feh_lst, alpha_lst)

fnew = feh_to_z(feh, alpha)
cax = ax1.imshow(fnew, interpolation='nearest', aspect='auto')
y1, y2 = ax1.get_ylim()
ax1.set_ylim(y2, y1)
ax1.set_xticklabels([v*dfe+fe0 for v in ax1.get_xticks()])
ax1.set_yticklabels([v*daf+af0 for v in ax1.get_yticks()])

ax2.plot_surface(feh, alpha, fnew, cmap='jet', rstride=1, cstride=1,
                 lw=0, antialiased=True)
ax2.view_init(elev=30, azim=225)

family = 'Georgia'
cbar = fig.colorbar(cax, cax=ax3)
cbar.set_label('Z', family=family)
for ax in [ax1, ax2]:
    ax.set_xlabel('[Fe/H]', family=family)
    ax.set_ylabel(u'[\u03b1/Fe]', family=family)
ax2.set_zlabel('Z', family=family)
ax2.set_xticks(np.arange(fe0, fe1+1e-6, 1.0))

fig.savefig('zmetal.png')
plt.show()
