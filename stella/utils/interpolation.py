def parabolic(ax,ay,x):
    '''
    parabolic interpolation
    if x = ax[1], f and h stays, but parab becomes
    parab = ay[0]+(ax[1]-ax[0])*f
          = ay[0]+(ay[1]-ay[0])
          = ay[1]
    '''
    f = (ay[1]-ay[0])/(ax[1]-ax[0])
    h = ((ay[2]-ay[0])/(ax[2]-ax[0])-f)/(ax[2]-ax[1])
    parab = ay[0]+(x-ax[0])*(f+h*(x-ax[1]))
    return parab

def newton(ax,ay,x):
    '''
    newton interpolation
    '''
    n = len(ax)
    y = [v for v in ay]

    for k in xrange(n-1):
        for l in xrange(n-k-1):
            y[l] = (y[l+1]-y[l])/(ax[l+k+1]-ax[l])

    ybar = y[0]
    for m in xrange(1,n):
        ybar = ybar*(x-ax[m])+y[m]
    return ybar
