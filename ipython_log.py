# IPython log file

import rings
R = rings.RingSetManager("LFI28M", 128, "full", fixfactor=1e3)
config=dict(input_map="/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits")
from planck import Planck
ch = Planck()["LFI28M"]
config["nside"]=128
get_ipython().magic('log')
get_ipython().magic('log start')
get_ipython().magic('logstart')
input_map = np.array(hp.ud_grade(
                                    hp.read_map(config["input_map"] % ch.f.freq, (0,1,2),  nest=True),
                                    config["nside"],
                                    order_in="NESTED",
                                    order_out="NESTED"
                                    )
                        ) * 1e3
import numpy as np
import healpy as hp
input_map = np.array(hp.ud_grade(
                                    hp.read_map(config["input_map"] % ch.f.freq, (0,1,2),  nest=True),
                                    config["nside"],
                                    order_in="NESTED",
                                    order_out="NESTED"
                                    )
                        ) * 1e3
pd.Series(input_map[0])
import pandas as pd
s = pd.Series(input_map[0])
s
s.index
s.index.name
R.data.pix
R.data.pix.max()
s.index.max()
s_rings = s.reindex(R.data.pix)
