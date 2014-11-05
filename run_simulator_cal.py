from planck import Planck
pl=Planck()
from testenv import cluster
#from planck import private

freqs = 30, 44
freqs = 70,
freqs = 30, 44, 70

def make_config_string(c):
    return " ".join(["--%s='%s'" % (k,str(v)) for k,v in c.items()])

config = dict( 
nside = 64,
#mask_filename = "largemask_%d.fits", # % ch.f.freq
mask_filename = "/global/project/projectdirs/planck/software/zonca/dev/chi2cal/destripingmask_%d.fits",
only_orb_dip = False,
destripe = True,
precond = True,
#tag = "3-%d" % private.survey[5].PID_LFI[1],
tag = "full",
ddx9data = False,
straylight = False,
pencil = False,
input_map = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits",
#input_map = "",
preremove_sol_dip=False,
remove_polarization=False,
init_dipole_fit=False,
scale_sol_dip_straylight=1.,
datarelease="dx11_delta",
#input_cal = "DX11D",
input_cal = "DX11DSLOW",
white_noise_scale = 1.,
input_map_polarization = False,
)

config_strings = []

#config_strings.append(make_config_string(config))
config["input_map_polarization"] = False
#config_strings.append(make_config_string(config))
config["scale_sol_dip_straylight"] = .8
config_strings.append(make_config_string(config))

for config_string in config_strings:
    for freq in freqs:
        for ch in pl.f[freq].ch:
            cluster.run_serial("itcal_%s" % ch.tag, "python simulator_cal.py --chtag {chtag} ".format(chtag=ch.tag) + config_string, mem=20)
