#!/usr/bin/env python3
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import SkyCoord
from stella.catalog  import find_catalog
from stella.kinetics import potential
from stella.kinetics import orbit
from stella.constant import pc


solar_uvw = (9.6, 255.2, 9.3) # from Reid et al., 2014
R0 = 8.34                     # from Reid et al., 2014

potential_lst = [potential.PointPotential(M=4.3e6), # from Gillessen et al. 2009
                 potential.HernquistPotential(M=3.76e9, a=0.1), # from Kenyon et al. 2014
                 potential.MiyamotoNagaiPotential(M=6e10, a=2.75, b=0.3),
                 potential.NFWPotential(M=1e12, rs=20),
                ]
vcirc_lst = list(map(lambda potential: potential.v_circ(R0), potential_lst))
v0 = np.sqrt((np.array(vcirc_lst)**2).sum())
print(vcirc_lst, v0)

t_lst = np.arange(0, 0.5, 0.0001) # in Gyr

x_lst, y_lst, z_lst = orbit.compute_Galorbit(
                        potential = potential_lst,
                        xyz=(R0,0.,0.),
                        uvw=(0.,0.,0.),
                        solar_uvw=solar_uvw,
                        t=t_lst)

# HD 122563 = HIP 68594
hip = 68594

item = find_catalog.find_HIP2(hip)
ra = item['RAdeg']
dec = item['DEdeg']
rv = (-26.58, 0.15)
parallax = (item['Plx'], item['e_Plx'])
pm = ((item['pmRA'], item['e_pmRA']),(item['pmDE'], item['e_pmDE']))
uvw = orbit.compute_UVW(ra=ra,dec=dec,rv=rv,parallax=parallax,pm=pm,U_plus='center')
xyz = orbit.compute_GalXYZ(ra=ra,dec=dec,parallax=parallax,R0=R0)
x1_lst, y1_lst, z1_lst = orbit.compute_Galorbit(
                        potential = potential_lst,
                        xyz=xyz,
                        uvw=uvw,
                        solar_uvw=solar_uvw,
                        t=t_lst)

def plot(t_lst, x_lst, y_lst, z_lst, color, filename):

    r_lst = np.sqrt(x_lst**2 + y_lst**2 + z_lst**2)

    fig1 = plt.figure(figsize=(16,4),dpi=150)
    ax1 = fig1.add_axes([0.04,0.12,0.19,0.8])
    ax2 = fig1.add_axes([0.28,0.12,0.19,0.8])
    ax3 = fig1.add_axes([0.52,0.12,0.19,0.8])
    ax4 = fig1.add_axes([0.76,0.12,0.2,0.8])

    ax1.plot(t_lst, x_lst, color+'-')
    ax1.plot(t_lst, y_lst, color+'--')
    ax2.plot(t_lst, z_lst, color+':')
    ax3.plot(t_lst, r_lst, color+'-')
    ax4.plot(x_lst, y_lst, color+'-')
    ax4.axhline(color='k',ls='--')
    ax4.axvline(color='k',ls='--')

    for ax in [ax1, ax2, ax3, ax4]:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(9)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(9)
        ax.set_xlabel('t (Gyr)',fontsize=9)
    ax1.set_ylabel('X, Y (kpc)',fontsize=9)
    ax2.set_ylabel('Z (kpc)',fontsize=9)
    ax3.set_ylabel('r (kpc)',fontsize=9)
    ax4.set_xlabel('X (kpc)',fontsize=9)
    ax4.set_ylabel('Y (kpc)',fontsize=9)
    
    fig1.savefig(filename)
    plt.close(fig1)


plot(t_lst, x_lst, y_lst, z_lst, 'r', 'orbit_sun.png')
plot(t_lst, x1_lst, y1_lst, z1_lst, 'b', 'orbit_HD122563.png')


fig = plt.figure(dpi=150)
ax = fig.gca(projection='3d')
ax.plot(x_lst, y_lst, z_lst, 'r-')
ax.plot(x1_lst, y1_lst, z1_lst, 'b-')
x1, x2 = ax.get_xlim()
y1, y2 = ax.get_ylim()
z1, z2 = ax.get_zlim()
ax.plot([x1,x2],[0,0],[0,0],'k-',alpha=0.5)
ax.plot([0,0],[y1,y2],[0,0],'k-',alpha=0.5)
ax.plot([0,0],[0,0],[z1,z2],'k-',alpha=0.5)
ax.set_xlim(x1,x2)
ax.set_ylim(y1,y2)
ax.set_zlim(z1,z2)
#ax.plot(x_lst, y_lst, 'r-',zs=z1,zdir='z',alpha=0.5)
ax.set_xlabel('X (kpc)')
ax.set_ylabel('Y (kpc)')
ax.set_zlabel('Z (kpc)')
plt.show()
