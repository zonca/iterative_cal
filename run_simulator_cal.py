from planck import Planck
pl=Planck()
from testenv import cluster
#from planck import private

freqs = 70,
freqs = 30, 44, 70
freqs = 44,
freqs = 70,
freqs = 30, 44, 70

def make_config_string(c):
    return " ".join(["--%s='%s'" % (k,str(v)) for k,v in c.items()])

config = dict( 
nside = 128,
#mask_filename = "largemask_%d.fits", # % ch.f.freq
only_orb_dip = False,
destripe = True,
precond = True,
#tag = "3-%d" % private.survey[5].PID_LFI[1],
tag = "full",
ddx9data = False,
straylight = False,
unknown_straylight = False,
pencil = True,
input_map = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits",
#input_map = "",
preremove_sol_dip=False,
remove_polarization=False,
init_dipole_fit=False,
scale_sol_dip_straylight=1.,
datarelease="dx11_delta",
input_cal = "DX11D",
#input_cal = "DX11DDVV",
#input_cal = "WAT1",
#input_cal = "DX11DSLOW",
white_noise_scale = 1.,
input_map_polarization = False,
dipole_constraint = "sol_dip",
remove_dipoles_signal = False
)

config_strings = []

#config_strings.append(make_config_string(config))
config["scale_sol_dip_straylight"] = 1.0
config["correct_main_beam_eff"] = False
config["mask_filename"] = "/global/project/projectdirs/planck/software/zonca/dev/chi2cal/destripingmask_30.fits"
#config["input_conv_sol_dip"] = True
#config["input_conv_orb_dip"] = True
config["straylight"] = False
config["scale_straylight"] = 1
#config_strings.append(make_config_string(config))
config["unknown_straylight"] = False

config["pencil"] = False
config["pencil_input"] = False
#config["dipole_constraint"] = ""
n = "S_DX11D_C"

config_strings.append((n, make_config_string(config)))

#for n in [1,2,3]:
for name,config_string in config_strings:
    for freq in freqs:
        for ch in pl.f[freq].ch:
            #print("python simulator_cal.py --chtag {chtag} ".format(chtag=ch.tag) + config_string)
            #if ch.tag == "LFI18M":
             cluster.run_serial(name + "_%s" % ch.tag, "python simulator_cal.py --chtag {chtag} ".format(chtag=ch.tag) + config_string, mem=20)
