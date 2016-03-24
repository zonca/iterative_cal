from scipy.sparse import linalg
import healpy as hp
import logging as l
import numpy as np
import datetime
from planck import Planck
from planck.pointing import compute_pol_weigths
import pandas as pd
import rings
import rings.calibration as rcal
import pymongo
#from planck import private
import argparse
import gc

parser = argparse.ArgumentParser(description='Iterative calibration')

l.basicConfig(level=l.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
l.root.level = l.DEBUG

config = dict( 
nside = 64,
mask_filename = "/global/project/projectdirs/planck/software/zonca/dev/chi2cal/destripingmask_%d.fits",
only_orb_dip = False,
destripe = True,
precond = True,
tag = "full",
ddx9data = False,
straylight = False,
unknown_straylight = False,
pencil = True,
input_map = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits",
preremove_sol_dip=False,
remove_polarization=False,
init_dipole_fit=False,
scale_sol_dip_straylight=1.,
datarelease="dx11_delta",
input_cal = "DX11DSLOW",
white_noise_scale = 1.,
input_map_polarization = False,
correct_main_beam_eff=False,
dipole_constraint = "sol_dip",
remove_dipoles_signal = False,
)

config = {'nside': 128, 'input_cal': 'DX11DDVV', 'destripe': True, 'init_dipole_fit': False, 'ddx9data': False, 'dipole_constraint': '', 'input_map_polarization': False, 'correct_main_beam_eff': False, 'pencil': False, 'remove_dipoles_signal': False, 'precond': True, 'datarelease': 'dx11_delta', 'unknown_straylight': False, 'preremove_sol_dip': False, 'only_orb_dip': False, 'remove_polarization': False, 'straylight': False, 
'input_map': '/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits', 'scale_sol_dip_straylight': 1.0, 'mask_filename': '/global/project/projectdirs/planck/software/zonca/dev/chi2cal/destripingmask_30.fits', 'tag': 'full', 'white_noise_scale': 0.0}
config["scale_straylight"] = 1.
config["chtag"] = "LFI28M"
config["remove_polarization_source"] = ""
config["pencil_input"] = True


fsl_fractional_solid_angle = {
"LFI18S": 0.62,
"LFI18M": 0.38,
"LFI19S": 0.58,
"LFI19M": 0.6 ,
"LFI20S": 0.7 ,
"LFI20M": 0.63,
"LFI21S": 0.7 ,
"LFI21M": 0.59,
"LFI22S": 0.5 ,
"LFI22M": 0.44,
"LFI23S": 0.43,
"LFI23M": 0.35,
"LFI24S": 0.15,
"LFI24M": 0.15,
"LFI25S": 0.06,
"LFI25M": 0.08,
"LFI26S": 0.05,
"LFI26M": 0.08,
"LFI27S": 0.76,
"LFI27M": 0.64,
"LFI28S": 0.83,
"LFI28M": 0.62,
}

def t_or_f(arg):
    ua = str(arg).upper()
    if ua == 'TRUE'[:len(ua)]:
       return True
    elif ua == 'FALSE'[:len(ua)]:
       return False
    else:
       pass  #error condition maybe?

for k,v in config.items():
    if type(v) is bool:
        argtype = t_or_f
    else:
        argtype = type(v)
    parser.add_argument("--" + k, type=argtype) 

args = vars(parser.parse_args())

for k,v in config.items():
    if not args[k] is None:
        config[k] = args[k]

config["versions"] = dict(rings="1.0", iterative_cal="1.0", planck="1.0")
#config["versions"] = dict(rings=get_module_gitcommit(rings), iterative_cal=get_folder_gitcommit("."), planck=get_module_gitcommit(planck))

ch = Planck()[config["chtag"]]

if not config["precond"]:
    prec = "np"
else:
    prec = ""

R = rings.RingSetManager(ch.tag, config["nside"], tag=config["tag"], by_ring=True, del_psi=False, ringsets_folder="/global/scratch2/sd/planck/user/zonca/data/ringsets_%s" % config["datarelease"], fixfactor=1e3)

if config["mask_filename"]:
    try:
        R.apply_mask(config["mask_filename"] % ch.f.freq)
    except:
        R.apply_mask(config["mask_filename"])

M = R.invert_invM(R.create_invM(R.data.index))

if config["dipole_constraint"]:
    print("Dipole contraint setup")
    if config["dipole_constraint"].startswith("conv"):
        _, dipole_map = R.destripe(R.data[config["dipole_constraint"]], maxiter=50, M=M)
    else:
        dipole_map = R.create_bin_map(R.data[config["dipole_constraint"]], M=M)
    dipole_map_cond = R.compute_dipole_constraint_invcond(M, dipole_map)
else:
    dipole_map = None
    dipole_map_cond = None

if config["input_map"] == "":
    R.data.c[:] = 0
elif not config["input_map"] == "data":
    print("Rescan input map")
    input_map = np.array(hp.ud_grade(
                                    hp.read_map(config["input_map"] % ch.f.freq, (0,1,2), nest=True),
                                    config["nside"],
                                    order_in="NESTED",
                                    order_out="NESTED"
                                    )
                        ) * 1e3
    R.data.c = pd.Series(input_map[0]).reindex(R.data.pix).values
    if config["input_map_polarization"]:
        qw, uw = compute_pol_weigths(R.data["psi"])
        R.data.c += pd.Series(input_map[1]).reindex(R.data.pix).values * qw
        R.data.c += pd.Series(input_map[2]).reindex(R.data.pix).values * uw

assert np.isnan(R.data.c).sum() == 0


if not config["input_map"] == "data":


    #if config["scale_sol_dip_straylight"] != 1.: 
    #    print("Scale sol dip straylight %.2f" % config["scale_sol_dip_straylight"])
    #    sol_dip_mainbeam = R.create_bin_map(R.data.sol_dip, M)

    #    sol_dip_straylight = R.remove_signal(R.data.sol_dip, bin_map=sol_dip_mainbeam)
    #    sol_dip_straylight -= sol_dip_straylight.mean()
    #    R.data.c += config["scale_sol_dip_straylight"] * sol_dip_straylight 

    #    if config["correct_main_beam_eff"]:
    #        main_beam_correction_factor = 1 + fsl_fractional_solid_angle[ch.tag]/100. * (1-config["scale_sol_dip_straylight"])
    #        sol_dip_mainbeam *= main_beam_correction_factor
    #        print("Scale main beam by %.4f" % (main_beam_correction_factor))
    #    R.data.c += sol_dip_mainbeam.I.reindex(R.data.pix).values
    #    del sol_dip_straylight
    #    del sol_dip_mainbeam
    #else:

    if config["pencil_input"]:
        conv = ""
    else:
        conv = "conv_"

    R.data.c += R.data[conv + "sol_dip"]
    R.data.c += R.data[conv + "orb_dip"]

    print("Galactic straylight")

    try:
        R.data["straylight"] *= config["scale_straylight"]
    except:
        pass

    if config["straylight"]: # add to data
        R.data.c += R.data.straylight
        if not config["unknown_straylight"]: # add to model
            R.data.orb_dip += R.data.straylight
    # 4pi input and pencil calibration, need to correct by beam efficiency
    if config["pencil"] and not config["pencil_input"]:
        beam_eff = pd.read_hdf("beam_efficiency_correction_factors.h5", "data")[ch.tag]
    else:
        beam_eff = 1.
    R.data.c *= beam_eff

    if config["white_noise_scale"]:
        white_noise_sigma = ch.white_noise_sigma * 1e6
        R.data.c += np.random.normal(size=len(R.data.hits)) * np.sqrt(white_noise_sigma * config["white_noise_scale"] / R.data.hits)

if not config["pencil"]:
    R.data.sol_dip = R.data.conv_sol_dip
    R.data.orb_dip = R.data.conv_orb_dip

try:
    del R.data["conv_sol_dip"]
    del R.data["conv_orb_dip"]
    del R.data["straylight"]
except:
    pass

print ("Decalibrate")
if config["input_cal"]:
    R.data.c /= rings.load_fits_gains(config["input_cal"], ch.tag, "DX10", by_ring=True).gain.reindex(R.data.index).fillna(method="ffill").fillna(method="bfill")

if config["only_orb_dip"]:
    R.data.sol_dip[:] = 0 #R.data.tot_dip - R.data.sol_dip - R.data.orb_dip 

if config["remove_dipoles_signal"]:
    R.data.sol_dip = R.remove_signal(R.data.sol_dip, M=M, dipole_map=dipole_map, dipole_map_cond=dipole_map_cond)
    R.data.orb_dip = R.remove_signal(R.data.orb_dip, M=M, dipole_map=dipole_map, dipole_map_cond=dipole_map_cond)

if config["remove_polarization"]:
    def remove_pol(R):
        qw, uw = compute_pol_weigths(R.data["psi"])
        pol_ringsets = pd.Series(0, index=R.data.index)
        pol_source = {
        "CMD": dict(
                folder = "/global/scratch2/sd/planck/user/peterm/commander_foregrounds/pol_default/",
                comps = ["synch", "dust"],
                to_mK = 1e-3,
                filename = "{comp}_{freq:03d}-{short_tag}_k00000.fits"),
        "CMDSHIFT": dict(
                folder = "/global/scratch2/sd/planck/user/zonca/commander_pol_bshift/",
                comps = ["synch", "dust"],
                to_mK = 1e-3,
                filename = "{comp}_{freq:03d}-{short_tag}_k00000.fits"),
        "DPC" : dict(
                folder = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/bandpass_corrected_maps_v2/",
                comps = [""],
                to_mK = 1e3,
                filename = "DX11D_{freq:03d}_128_conv120_corrected_full_v2.fits")
        }[config["remove_polarization_source"]]

        for comp in pol_source["comps"]:
            polarization_map = np.array(hp.ud_grade(hp.read_map(pol_source["folder"] + pol_source["filename"].format(comp=comp, freq=ch.f.freq, short_tag=ch.tag[3:], tag=ch.tag), (0,1,2), nest=True), config["nside"], order_in="NEST", order_out="NEST")) * pol_source["to_mK"]
            for i,w in [(1, qw), (2, uw)]:
                pol_ringsets += pd.Series(polarization_map[i]).reindex(R.data.pix).values * w

        pol_ringsets /= rings.load_fits_gains("DX11D", ch.tag, "DX10", by_ring=True).gain.reindex(R.data.index).fillna(method="ffill").fillna(method="bfill")
        R.data.c -= pol_ringsets
    remove_pol(R)
    gc.collect()

del R.data["psi"]

if config["init_dipole_fit"]:
    out = pd.DataFrame(rings.calibration.dipole_fit(R, R.data.sol_dip+R.data.orb_dip), index=R.pids)
else:
    out = pd.DataFrame(dict(g0=1., o0=0.), index=R.pids)

initial_guess = np.concatenate([np.array(out["o0"]), np.array(out["g0"])])

num_outer_iterations = 1+6
n_ods = len(R.pids)

client = pymongo.MongoClient('mongodb01.nersc.gov')
db = client["planckcal"]
db.authenticate("planckcal", "dithoTX13")
gainsdb = db.gains

config["pids"] = R.pids.astype(np.float).tolist()

def store_results(iteration, gains, offsets, weights):
    if "_id" in config:
        del config["_id"]
    config.update(dict(
        iteration = iteration, 
        gains     = gains.tolist(), 
        offsets   = offsets.tolist(),
        weights   = weights.tolist()
        ))
    config["date"] = datetime.datetime.now()
    gainsdb.insert(config)

store_results(0, out.g0, out.o0, np.ones_like(out.g0))
g_prev = out.g0

hits_per_pp = rings.sum_by(R.data.hits, R.data.index, target_index=R.pids)

for outer_i in range(1, num_outer_iterations):
    print(("Iteration:", outer_i))
    #g_prev = pd.rolling_mean(outgains[outer_i-1], 5).fillna(method="bfill")

    if config["destripe"]:
        m_prev_bin, m_prev, baselines = R.destripe(R.data.c / g_prev.reindex(R.data.index) - R.data.orb_dip - R.data.sol_dip, return_baselines=True, maxiter=50, M=M)
        R.data.c[:] = R.remove_baselines(R.data.c, baselines * g_prev.reindex(baselines.index))
    else:
        m_prev_bin = R.create_bin_map(R.data.c / g_prev.reindex(R.data.index) - R.data.orb_dip - R.data.sol_dip, M=M)
        m_prev = m_prev_bin

    if config["dipole_constraint"]:
        m_prev -= R.fit_mono_dipole(m_prev, M, dipole_map, dipole_map_cond)
    RHS = rcal.create_RHS(R, g_prev, m_prev, M=M, dipole_map=dipole_map, dipole_map_cond=dipole_map_cond)
    matvec = rcal.create_matvec(R, g_prev, m_prev, M=M, dipole_map=dipole_map, dipole_map_cond=dipole_map_cond)
    LinCalOperator = linalg.LinearOperator(shape=(2*n_ods, 2*n_ods), matvec=matvec, dtype=np.double)
    Dinv = rcal.create_Dinv(R, hits_per_pp, m_prev)
    Dinv = rcal.mult_det(Dinv)
    if config["precond"]:
        precond_matvec = rcal.create_precon_matvec(Dinv)
    else:
        precond_matvec = rcal.create_preconv_matvec_hits(hits_per_pp)
    monitor = rcal.Monitor()
    LinCalPrecondOperator = linalg.LinearOperator(shape=(2*n_ods, 2*n_ods), matvec=precond_matvec, dtype=np.double)
    solution, info = linalg.cg(LinCalOperator, RHS, x0=initial_guess, tol=1e-15, maxiter=50, callback=monitor, M=LinCalPrecondOperator)

    g_prev = pd.Series(solution[n_ods:], index=R.pids)
    store_results(outer_i, g_prev, solution[:n_ods], 1./Dinv["11"])
    client.close()
