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
from planck import private
import pymongo
import argparse

parser = argparse.ArgumentParser(description='Iterative calibration')

l.basicConfig(level=l.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
l.root.level = l.DEBUG

config = dict( 
nside = 64,
chtag = "LFI19M",
mask_filename = "largemask_%d.fits", # % ch.f.freq
only_orb_dip = True,
destripe = False,
precond = True,
tag = "3-%d" % private.survey[5].PID_LFI[1],
ddx9data = False,
straylight = False,
pencil = False,
input_map = "",
datarelease = "dx10",
preremove_sol_dip = False,
remove_polarization = False,
init_dipole_fit = False,
scale_sol_dip_straylight=1.,
input_cal = "DX11D",
white_noise_scale = 0.,
input_map_polarization = False,
dipole_constraint = False,
)

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

R = rings.RingSetManager([ch.tag], config["nside"], tag=config["tag"], by_ring=True, del_psi=False, ringsets_folder="/global/scratch2/sd/planck/user/zonca/data/ringsets_%s" % config["datarelease"])

if config["mask_filename"]:
    R.apply_mask(config["mask_filename"] % ch.f.freq)

if config["input_map"] == "":
    R.data.c[:] = 0
else:
    input_map = np.array(hp.read_map(config["input_map"] % ch.f.freq, (0,1,2), nest=True)) * 1e3
    R.data.c = pd.Series(input_map[0]).reindex(R.data.index, level="pix")
    if config["input_map_polarization"]:
        qw, uw = compute_pol_weigths(R.data["psi"])
        R.data.c += pd.Series(input_map[1]).reindex(R.data.index, level="pix") * qw
        R.data.c += pd.Series(input_map[2]).reindex(R.data.index, level="pix") * uw

if config["white_noise_scale"]:
    white_noise_sigma = ch.white_noise_sigma * 1e6
    R.data.c += np.random.normal(size=len(R.data.hits)) * np.sqrt(white_noise_sigma * config["white_noise_scale"] / R.data.hits)

assert np.isnan(R.data.c).sum() == 0

M = R.invert_invM(R.create_invM(R.data.index))

hits_per_pp_series = R.data.hits.groupby(level="od").sum()
hits_per_pp = np.array(hits_per_pp_series)
R.pids = hits_per_pp_series.index
del hits_per_pp_series

if config["preremove_sol_dip"]:
    signal_removed_sol_dip = R.remove_signal(R.data.conv_sol_dip, M=M)
    signal_removed_sol_dip /= rings.load_fits_gains("PSEU2", ch.tag, "DX10", by_ring=True).gain.reindex(signal_removed_sol_dip.index, level="od").fillna(method="ffill").fillna(method="bfill")
    if config["scale_sol_dip_straylight"]:
        signal_removed_sol_dip *= config["scale_sol_dip_straylight"]
    R.data.c -= signal_removed_sol_dip
    del signal_removed_sol_dip

if config["dipole_constraint"]:
    if config["dipole_constraint"].startswith("conv"):
        _, dipole_map = R.destripe(R.data[config["dipole_constraint"]], maxiter=50, M=M)
    else:
        dipole_map = R.create_bin_map(R.data[config["dipole_constraint"]], M=M)
    dipole_map_cond = R.compute_dipole_constraint_invcond(M, dipole_map)

if not config["pencil"]:
    R.data.sol_dip[:] = R.data.conv_sol_dip
    R.data.orb_dip[:] = R.data.conv_orb_dip
del R.data["conv_sol_dip"]
del R.data["conv_orb_dip"]

R.data.c += R.data.sol_dip
if config["scale_sol_dip_straylight"]:
    R.data.c -= (1. - config["scale_sol_dip_straylight"]) * R.remove_signal(R.data.sol_dip)
R.data.c += R.data.orb_dip

print ("Decalibrate")
R.data.c /= rings.load_fits_gains(config["input_cal"], ch.tag, "DX10", by_ring=True).gain.reindex(R.data.index, level="od").fillna(method="ffill").fillna(method="bfill")

if config["straylight"]:
    R.data.orb_dip += R.data.straylight

#cmb = pd.Series(1e3*hp.ud_grade(hp.read_map("cmb_only_%dGHz.fits" % ch.f.freq), nside))
#R.data.c += cmb.reindex(R.data.index, level="pix")

if config["only_orb_dip"]:
    R.data.sol_dip[:] = 0 #R.data.tot_dip - R.data.sol_dip - R.data.orb_dip 

R.data.sol_dip = R.remove_signal(R.data.sol_dip, M=M)
R.data.orb_dip = R.remove_signal(R.data.orb_dip, M=M)

if config["remove_polarization"]:
    print ("Remove polarization")
    qw, uw = compute_pol_weigths(R.data["psi"])
    folder = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/lfi/BandPass_Corrected_Maps/"
    if not config["nside"] == 256:
        raise TypeError("remove_polarization works only at nside 256")
    polarization_map = np.array(hp.read_map(folder + "LFI_SkyMap-bandpassCorrected_%03d_0256_DX10_full.fits" % ch.f.freq, (0,1,2), nest=True)) * 1e3
    pol_ringsets =  pd.Series(polarization_map[1]).reindex(R.data.index, level="pix") * qw
    pol_ringsets += pd.Series(polarization_map[2]).reindex(R.data.index, level="pix") * uw
    pol_ringsets /= rings.load_fits_gains("PSEU2", ch.tag, "DX10", by_ring=True).gain.reindex(pol_ringsets.index, level="od").fillna(method="ffill").fillna(method="bfill")
    R.data.c -= pol_ringsets
    del polarization_map, qw, uw, pol_ringsets

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


for outer_i in range(1, num_outer_iterations):
    print(("Iteration:", outer_i))
    #g_prev = pd.rolling_mean(outgains[outer_i-1], 5).fillna(method="bfill")

    if config["destripe"]:
        m_prev_bin, m_prev, baselines = R.destripe(R.data.c / g_prev.reindex(R.data.index, level="od") - R.data.orb_dip - R.data.sol_dip, return_baselines=True, maxiter=50, M=M)
        R.data.c[:] = R.remove_baselines(R.data.c, baselines * g_prev.reindex(baselines.index, level="od"))
    else:
        m_prev_bin = R.create_bin_map(R.data.c / g_prev.reindex(R.data.index, level="od") - R.data.orb_dip - R.data.sol_dip, M=M)
        m_prev = m_prev_bin
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
