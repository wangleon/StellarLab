#!/usr/bin/env python3
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from astropy.coordinates import SkyCoord
from stella.catalog  import HIP2
from stella.kinetics import potential
from stella.kinetics import orbit
from stella.constant import pc

def main():

    solar_uvw = (9.6, 255.2, 9.3) # from Reid et al., 2014
    R0 = 8.34                     # from Reid et al., 2014
    
    potential_lst = [potential.PointPotential(M=4.3e6), # from Gillessen et al. 2009
                     potential.HernquistPotential(M=3.76e9, a=0.1), # from Kenyon et al. 2014
                     potential.MiyamotoNagaiPotential(M=6e10, a=2.75, b=0.3),
                     potential.NFWPotential(M=1e12, rs=20),
                    ]
    vcirc_lst = list(map(lambda potential: potential.get_vcirc(R0), potential_lst))
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
    item = HIP2.find_object(hip)
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
    
    
    fig = plt.figure(dpi=150)
    ax = fig.gca(projection='3d')
    #ax.plot(x_lst, y_lst, z_lst, 'r-')

    t_lst = t_lst[::10]
    x_lst = x_lst[::10]
    y_lst = y_lst[::10]
    z_lst = z_lst[::10]
    x1_lst = x1_lst[::10]
    y1_lst = y1_lst[::10]
    z1_lst = z1_lst[::10]

    x1, x2 = -10, 10
    y1, y2 = -10, 10
    z1, z2 = -1, 1
    def update(i):
        ax.cla()
        ax.set_xlim(x1, x2)
        ax.set_ylim(y1, y2)
        ax.set_zlim(z1, z2)
        ax.set_title('t=%5.2f Gyr'%t_lst[i])

        ax.plot([x1,x2],[0,0],[0,0],'k-',alpha=0.5, lw=0.5)
        ax.plot([0,0],[y1,y2],[0,0],'k-',alpha=0.5, lw=0.5)
        ax.plot([0,0],[0,0],[z1,z2],'k-',alpha=0.5, lw=0.5)
        line, = ax.plot(x_lst[0:i], y_lst[0:i], z_lst[0:i], 'r-')
        line, = ax.plot(x1_lst[0:i], y1_lst[0:i], z1_lst[0:i], 'b-')
        ax.view_init(elev=30, azim=360*3/t_lst.size-i)
        ax.set_xlabel('X (kpc)')
        ax.set_ylabel('Y (kpc)')
        ax.set_zlabel('Z (kpc)')
        ax.grid(False)
        return line,

    anim = animation.FuncAnimation(fig, update, 
            frames=np.arange(t_lst.size),interval=1, blit=True)
    anim.save('orbits.mp4', fps=25, extra_args=['-vcodec', 'libx264'])

    #plt.show()

if __name__=='__main__':
    main()
