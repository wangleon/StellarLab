
class ColorIndexError(Exception):
    def __init__(self, index, ref):
        self.index = index
        self.ref   = ref
    def __str__(self):
        return 'Color index %s cannot be used for Teff calibration in %s'%(
                self.index, self.ref)

class ParamRangeError(Exception):
    def __init__(self, param, value, ref):
        self.param = param
        self.value = value
        self.ref   = ref
    def __str__(self):
        return 'Parameter (%s) = %g is not in the range of %s'%(
                self.param, self.value, self.ref)

class ApplicableRangeError(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return 'Input parameters excess the applicable range.'

class ParamMissingError(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return 'Missing Input parameters.'

class MissingParamError(Exception):
    def __init__(self, param, ref):
        self.param = param
        self.ref   = ref
    def __str__(self):
        return 'Missing input parameter %s in %s'%(self.param, self.ref)
