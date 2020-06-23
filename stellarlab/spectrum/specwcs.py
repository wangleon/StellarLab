import numpy as np
import numpy.polynomial as poly
import astropy.io.fits as fits

class MultiSpecItem(object):
    '''Class for MultiSpec format
    
    '''

    def __init__(self,string):
        self.string = string
        self.parse_string()

    def parse_string(self):
        # fix "xxxE-40.xxx problem"

        g = self.string.split()
        for i in range(len(g)):
            if g[i].count('.')==2:
                ## find first dot
                #pos1 = e.find(".")
                ## find second dot
                #pos2 = e[pos1+1:].find(".") + pos1 + 1
                ## replace the second dot
                #g[i] = g[i].replace('0.',' 0.')

                t = g[i].split('.')
                t[2] = t[1][-1]+'.'+t[2]
                t[1] = t[1][:-1]
                g[i] = t[0]+'.'+t[1]+' '+t[2]

        string = ' '.join(g)

        g = string.split()
        self.ap    =   int(g[0])
        self.beam  =   int(g[1])
        self.dtype =   int(g[2])
        self.w1    = float(g[3])
        self.dw    = float(g[4])
        self.nw    =   int(g[5])
        self.z     = float(g[6])
        self.aplow = float(g[7])
        self.aphigh= float(g[8])

        if self.dtype == 0:
            # linear coordinate
            pass

        elif self.dtype == 2:
            # nonlinear coordiante
            func = 0
            self.wt    = {}
            self.wt0   = {}
            self.ftype = {}
            self.order = {}
            self.pmin  = {}
            self.pmax  = {}
            self.c     = {}
            pos = 0
            while(pos+9<len(g)):
                self.wt[func]    = float(g[ 9+pos])
                self.wt0[func]   = float(g[10+pos])
                self.ftype[func] = round(float(g[11+pos]))
                if self.ftype[func] in [1,2]:
                    # Chebyshev (ftype=1) or Legendre (ftype=2) Polynomial
                    self.order[func] = int(g[12+pos])
                    self.pmin[func]  = float(g[13+pos])
                    self.pmax[func]  = float(g[14+pos])
                    self.c[func]     = [float(g[15+pos+j])
                                        for j in range(self.order[func])]
                    pos += 3+3+self.order[func]
                func += 1

    @staticmethod
    def get_wat2(head):
        '''Convert WAT2 string to a dict containing :class:`MultiSpecItem`
        instaces.

        Args:
            head (string): WAT2 string.
        Returns:
            dict: A dict containing :class:`MultiSpecItem` instances.
        '''
        string = ''
        for item in head.items():
            if item[0][0:4]=='WAT2':
                string += item[1].ljust(68)
        ap_num = head['NAXIS2']
        g = string.split('spec')
        pos_list = []
        for ap in range(1, ap_num+1):
            label = 'spec%d = '%ap
            pos_list.append(string.find(label))
        pos_list.append(-1)

        multispec_items = {}
        for ap in range(1, ap_num+1):
            pos1 = pos_list[ap-1]
            pos2 = pos_list[ap]
            tmp  = string[pos1:pos2].split('=')[1].replace('"','').strip()
            multispec_items[ap] = MultiSpecItem(tmp)
        return multispec_items

    def get_wv(self):
        '''Get wavelength for a record in multispec.

        Returns:
            :class:`numpy.array`: An array of wavelengths
        '''
        if self.dtype == 0:
            # linear coordinate
            wv = self.w1 + np.arange(self.nw)*self.dw

        elif self.dtype == 2:
            # nonlinear coordiante
            p = np.arange(self.pmin[0], self.pmax[0]+1e-6)
            n = (p - (self.pmax[0] + self.pmin[0])/2.0)/(
                     (self.pmax[0] - self.pmin[0])/2.0)
            wv  = np.zeros_like(n)
            for func in self.ftype.keys():
                if int(self.ftype[func])==1:
                    # ftype = 1: Chebyshev Polynomial
                    wvi = poly.chebyshev.chebval(n,self.c[func])
                elif int(self.ftype[func])==2:
                    # ftype = 2: Legendre Polynomial
                    wvi = poly.legendre.legval(n,self.c[func])
                elif int(self.ftype[func])==3:
                    # ftype = 3: Cubic Spline
                    pass
                wv += self.wt[func]*wvi
        return wv

def load_multispec(filename, key='beam'):
    '''Read the *Multispec* spectral world coordinate system.

    The dispersion functions are specified by strings with identifiers
    *"specN"*, where *N* is the phyiscal image line.
    The strings contain a series of fields::

        specN = ap beam dtype w1 dw nw z aplow aphigh [functions_i]

    where

        * **ap** - the aperture number.
        * **beam** - the beam number.
        * **dtype** - the dispersion type with the following meanings:

            ===== ===================================================================================
            Value Meaning
            ===== ===================================================================================
            -1    coordinates are not dispersion coordinates or spectrum is not dispersion calibrated
            0     linear dispersion sampling
            1     log-linear dispersion sampling
            2     nonlinear dispersion
            ===== ===================================================================================

        * **w1** - the dispersion coordinate of the first physical pixel.
        * **dw** - the average dispersion interval per physical pixel.
        * **nw** - the number of valid pixels.
        * **z** - the Doppler factor applied to all dispersion coordinates by
          multiplying by 1/(1 + *z*). 0 means no Doppler correction.
        * **aplow** - the lower limit of the extracted apertures.
        * **aphigh** - the upper limit of the extracted apertures.
        * **[functions_i]** - zero or more function descriptions. There are no
          such descriptions for linear or log-linear dispersion coordinate
          systems.
          For the nonlinear dispersion systems, the descriptions have the
          following fileds::

              function_i = wt_i w0_i ftype_i [parameters] [coefficients]

          where

          * **wt_i** - weight.
          * **w0_i** - zero point offset.
          * **ftype_i** - the function type codes with the following meanings:

            ===== =========================
            Value Meaning
            ===== =========================
            1     Chebyshev polynomial
            2     Legendre polynomial
            3     Cubic spline
            4     Linear spline
            5     Pixel coordinate array
            6     Sampled coordinate array
            ===== =========================

          * **[parameters]** - the parameters of the dispersion function 
          * **[coefficients]** - the coefficients of the dispersion function

          The final wavelength is the weighted sum of *nfunc* individual
          dispersion functions *W*:sub:`i`\ (*p*):

          .. math::

              w=\sum_{i=1}^{\mathrm{nfunc}}\left(\\frac{w_{t,i}(w_{0,i}+W_i(p))}{1+z}\\right)


    Args:
        filename (string): The name of file to read.
        key (string): Key of the returned dict. Either "beam" or "ap".
    Returns:
        dict: A dict containing wavelength and flux of each order as tuples.
    Examples:

        .. code:: python

           import matplotlib.pyplot as plot
           from stella.spectrum import specwcs
           res = specwcs.load_multispec('HD2796.fits')
           # 'HD2796.fits' is a multi_spec fits

           fig = plt.figure()
           ax = fig.gca()
           for order, (wave, flux) in sorted(res.items()):
               ax.plot(wave, flux, '-')
           plt.show()
    '''

    data, head = fits.getdata(filename,header=True)
    multispec_items = MultiSpecItem.get_wat2(head)

    result = {}

    for N in multispec_items.keys():
        ap   = multispec_items[N].ap
        beam = multispec_items[N].beam
        wave = multispec_items[N].get_wv()
        if multispec_items[N].dtype==0:
            flux = data[N-1,:]
        elif multispec_items[N].dtype==2:
            pmin = multispec_items[N].pmin[0]
            pmax = multispec_items[N].pmax[0]
            flux = data[N-1,int(pmin)-1:int(pmax)]

        if key=='beam':
            result[beam] = (wave, flux)
        elif key=='ap':
            result[ap] = (wave, flux)
        else:
            print('Unknown key')
            raise ValueError
    return result
