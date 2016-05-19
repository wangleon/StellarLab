import numpy as np

def feh_to_z(feh,alpha=0.0):
    '''
    Convert stellar [Fe/H] and [alpha/Fe] to Z.
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
    return griddata(coor,val,(feh,alpha),method='cubic')
