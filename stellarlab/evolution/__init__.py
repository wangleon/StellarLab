from . import base
from .y2     import Y2
#from .y2pms  import Y2PMSTrack
#from .geneva import GenevaTrack
from . import geneva
from .yapsi import YaPSI

def get_track(track,**kwargs):
    track = track.lower().strip()
    if track == 'y2':
        return Y2Track(**kwargs)
    elif track == 'y2pms':
        return Y2PMSTrack(**kwargs)
    elif track == 'geneva':
        return GenevaTrack(**kwargs)
    else:
        raise ValueError
