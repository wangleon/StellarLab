import numpy as np

def feh_to_z(FeH, alpha=0.0):
    '''
    Convert stellar [Fe/H] abundances between [−3.0, +1.0] and [α/Fe] abundances
    between [0.0, +0.6] to *Z*.

    Args:
        FeH (float): [Fe/H] abundance ratio
        alpha (float, optional): [α/Fe] abundance ratio. Defaut is 0.0
    Returns:
        float: Metal ratio (*Z*)

    The returned values are shown as below:

    .. plot::

       import numpy as np
       import matplotlib.pyplot as plt
       from mpl_toolkits.mplot3d import Axes3D
       from stella.parameter.metal import feh_to_z

       fig = plt.figure(figsize=(10,4), dpi=200)
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

       plt.show()

    '''
    from scipy.interpolate import griddata
    ydata,xdata = np.mgrid[0.0:0.6+1e-3:0.3,-3.0:1.0+1e-3:0.5]
    coor = np.concatenate((xdata.reshape(-1,1),
                           ydata.reshape(-1,1)),axis=1)
    #z nodes
    z0 = np.array([
            [0.000019, 0.000062, 0.000195,
             0.000615, 0.001935, 0.006021,
             0.018120, 0.049711, 0.110798],
            [0.000032, 0.000102, 0.000321,
             0.001012, 0.003174, 0.009774,
             0.028557, 0.072793, 0.142689],
            [0.000058, 0.000182, 0.000574,
             0.001807, 0.005627, 0.016990,
             0.047000, 0.106471, 0.177489]
            ])
    val = z0.reshape(-1)
    return griddata(coor,val,(FeH,alpha),method='cubic')



def get_feh(**kwargs):
    ref = kwargs.pop('ref').strip().lower()

    if ref in ['onehag2009']:
        spectype = kwargs.pop('spectype')
        if spectype[0] == 'G':
            FeH = __get_gstar_feh_Onehag2009(**kwargs)
        elif spectype[0] == 'F':
            FeH = __get_fstar_feh_Onehag2009(**kwargs)

def __get_gstar_feh_Onehag2009(**kwargs):
    """metallicity calibration for G-stars
    based on calibration of Onehag, et al. 2009, A&A, 498, 527
    (2009A&A...498..527O) Appdendix A.3
    only applicable for logg >= 3.0
    """

    b_y = kwargs.pop('b-y')
    m1  = kwargs.pop('m1')
    c1  = kwargs.pop('c1')

    if 0.37 <=b_y <= 0.59 and \
       0.03 <= m1 <= 0.57 and \
       0.1  <= c1 <= 0.47:

        FeH = -2.796 + 39.21*m1 - 88.97*m1**2 - 73.43*m1*b_y \
              +181.4*m1**2*b_y + (27.03*m1 - 1.220*c1 - 41.42*m1**2)*c1
        if -2.6 <= FeH <= 0.4:
            return FeH
        else:
            raise ValueError

def __get_fstar_feh_Onehag2009(**kwargs):
    """metallicity calibration for F-stars
    based on calibration of Onehag, et al. 2009, A&A, 498, 527
    (2009A&A...498..527O) Appdendix A.4
    """

    b_y = kwargs.pop('b-y')
    m1  = kwargs.pop('m1')
    c1  = kwargs.pop('c1')

    if 0.22 <=b_y <= 0.38 and \
       0.03 <= m1 <= 0.21 and \
       0.17 <= c1 <= 0.58:
        c3 = 0.4462 - 2.233*b_y + 2.885*b_y**2
        FeH = 1.850 - 34.21*m1 + 105.43*m1*b_y \
              + 179.8*m1**2*b_y - 242.4*m1*b_y**2 \
              +(2.757-20.38*m1+0.2777*b_y)*math.log10(m1-c3)
        if -3.5 <= FeH <= 0.2:
            return FeH
        else:
            raise ValueError

