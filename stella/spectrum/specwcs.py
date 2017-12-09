def load_multispec(filename):
    wv_group, flux_group = specwcs.calib(filename)

    spec = Spec()
    for ord in wv_group.keys():
        specdata = SpecData(wv   = wv_group[ord],
                            flux = flux_group[ord])
        spec.add_order(ord,specdata)

    return spec


class specwcs(object):

    @staticmethod
    def calib(filename):

        data,head = pf.getdata(filename,header=True)
        multispec_items = MultiSpecItem.get_wat2(head)

        wv_group   = {}
        flux_group = {}

        for N in multispec_items.keys():
            ord = multispec_items[N].beam
            wv = multispec_items[N].get_wv()
            if multispec_items[N].dtype==0:
                flux = data[N-1,:]
            elif multispec_items[N].dtype==2:
                pmin = multispec_items[N].pmin[0]
                pmax = multispec_items[N].pmax[0]
                flux = data[N-1,int(pmin)-1:int(pmax)-int(pmin)+1]

            wv_group[ord]   = wv
            flux_group[ord] = flux

        return wv_group,flux_group

class MultiSpecItem(object):
    """Object for MultiSpec Format"""

    def __init__(self,string):
        self.string = string
        self.parse_string()

    def parse_string(self):
        # fix "xxxE-40.xxx problem"

        g = self.string.split(" ")
        for i in np.arange(len(g)):
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
            pos       = 0
            while(pos+9<len(g)):
                self.wt[func]    = float(g[ 9+pos])
                self.wt0[func]   = float(g[10+pos])
                self.ftype[func] = round(float(g[11+pos]))
                if self.ftype[func] in [1,2]:
                    # Chebyshev (ftype=1) or Legendre (ftype=2) Polynomial
                    self.order[func] = int(g[12+pos])
                    self.pmin[func]  = float(g[13+pos])
                    self.pmax[func]  = float(g[14+pos])
                    self.c[func]     = []
                    for j in np.arange(self.order[func]):
                        self.c[func].append(float(g[15+pos+j]))
                    pos += 3+3+self.order[func]
                func += 1

    @staticmethod
    def get_wat2(head):
        ''' get WAT2 string, return a turple, consist of
        MultiSpecItem instances '''
        string = ""
        for item in head.items():
            if item[0][0:4]=="WAT2":
                string += item[1].ljust(68)
        ap_num = head['NAXIS2']
        g = string.split('spec')
        pos_list = []
        for ap in np.arange(1,ap_num+1):
            label = "spec"+str(ap)+" = "
            pos_list.append(string.find(label))
        pos_list.append(-1)

        multispec_items = {}
        for ap in np.arange(1,ap_num+1):
            pos1 = pos_list[ap-1]
            pos2 = pos_list[ap]
            tmp  = string[pos1:pos2].split("=")[1].replace("\"","").strip()
            multispec_items[ap]=MultiSpecItem(tmp)
        return multispec_items

    def get_wv(self):
        ''' get wavelength for a record in multispec'''
        if self.dtype == 0:
            # linear coordinate
            wv = self.w1 + np.arange(self.nw)*self.dw

        elif self.dtype == 2:
            # nonlinear coordiante
            p   = np.arange(self.pmin[0], self.pmax[0]+1e-6)
            n   = (p - (self.pmax[0] + self.pmin[0])/2.0)/((
                        self.pmax[0] - self.pmin[0])/2.0)
            wv  = np.zeros_like(n)
            for func in self.ftype.keys():
                if int(self.ftype[func])==1:
                    # ftype = 1: Chebyshev Polynomial
                    wvi = chebyshev_poly(n,self.c[func])
                elif int(self.ftype[func])==2:
                    # ftype = 2: Legendre Polynomial
                    wvi = legendre_poly(n,self.c[func])
                elif int(self.ftype[func])==3:
                    # ftype = 3: Cubic Spline
                    pass
                wv += self.wt[func]*wvi
        return wv
    @staticmethod
    def calib(filename):

        data,head = pf.getdata(filename,header=True)
        multispec_items = MultiSpecItem.get_wat2(head)

        wv_group   = {}
        flux_group = {}

        for N in multispec_items.keys():
            ord = multispec_items[N].beam
            wv = multispec_items[N].get_wv()
            if multispec_items[N].dtype==0:
                flux = data[N-1,:]
            elif multispec_items[N].dtype==2:
                pmin = multispec_items[N].pmin[0]
                pmax = multispec_items[N].pmax[0]
                flux = data[N-1,int(pmin)-1:int(pmax)-int(pmin)+1]

            wv_group[ord]   = wv
            flux_group[ord] = flux

        return wv_group,flux_group

class MultiSpecItem(object):
    """Object for MultiSpec Format"""

    def __init__(self,string):
        self.string = string
        self.parse_string()

    def parse_string(self):
        # fix "xxxE-40.xxx problem"

        g = self.string.split(" ")
        for i in np.arange(len(g)):
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
            pos       = 0
            while(pos+9<len(g)):
                self.wt[func]    = float(g[ 9+pos])
                self.wt0[func]   = float(g[10+pos])
                self.ftype[func] = round(float(g[11+pos]))
                if self.ftype[func] in [1,2]:
                    # Chebyshev (ftype=1) or Legendre (ftype=2) Polynomial
                    self.order[func] = int(g[12+pos])
                    self.pmin[func]  = float(g[13+pos])
                    self.pmax[func]  = float(g[14+pos])
                    self.c[func]     = []
                    for j in np.arange(self.order[func]):
                        self.c[func].append(float(g[15+pos+j]))
                    pos += 3+3+self.order[func]
                func += 1

    @staticmethod
    def get_wat2(head):
        ''' get WAT2 string, return a turple, consist of
        MultiSpecItem instances '''
        string = ""
        for item in head.items():
            if item[0][0:4]=="WAT2":
                string += item[1].ljust(68)
        ap_num = head['NAXIS2']
        g = string.split('spec')
        pos_list = []
        for ap in np.arange(1,ap_num+1):
            label = "spec"+str(ap)+" = "
            pos_list.append(string.find(label))
        pos_list.append(-1)

        multispec_items = {}
        for ap in np.arange(1,ap_num+1):
            pos1 = pos_list[ap-1]
            pos2 = pos_list[ap]
            tmp  = string[pos1:pos2].split("=")[1].replace("\"","").strip()
            multispec_items[ap]=MultiSpecItem(tmp)
        return multispec_items

    def get_wv(self):
        ''' get wavelength for a record in multispec'''
        if self.dtype == 0:
            # linear coordinate
            wv = self.w1 + np.arange(self.nw)*self.dw

        elif self.dtype == 2:
            # nonlinear coordiante
            p   = np.arange(self.pmin[0], self.pmax[0]+1e-6)
            n   = (p - (self.pmax[0] + self.pmin[0])/2.0)/((
                        self.pmax[0] - self.pmin[0])/2.0)
            wv  = np.zeros_like(n)
            for func in self.ftype.keys():
                if int(self.ftype[func])==1:
                    # ftype = 1: Chebyshev Polynomial
                    wvi = chebyshev_poly(n,self.c[func])
                elif int(self.ftype[func])==2:
                    # ftype = 2: Legendre Polynomial
                    wvi = legendre_poly(n,self.c[func])
                elif int(self.ftype[func])==3:
                    # ftype = 3: Cubic Spline
                    pass
                wv += self.wt[func]*wvi
        return wv

