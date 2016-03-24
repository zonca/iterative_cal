import sqlite3
from planck.pointing import compute_pol_weigths
from scipy import interpolate
import rings
import pandas as pd
import pyfits
from planck import Planck
from planck.metadata import latest_exchange
from planck import private, pointing
import os
import numpy as np
import healpy as hp

def do1od(od, freq):
#od = 100
#freq = 30

    pl = Planck()

    NSIDE = 1024

    def get_4pi_dipole(obt, pix, psi, ch, orbtype):
        from dipole import Dipole
        dip = Dipole(obt=obt, type=orbtype) 
        theta, phi = hp.pix2ang(NSIDE, pix, nest=False)
        return dip.get_4piconv_dx10(ch, theta, phi, psi)

    print("-----OD %d" % od)

    # CONFIGURATION

    conf = {}
    #cal_tag = "SDX11DSLOW08SC1"
    #cal_tag = "SDX11DSLOWC1"
    #cal_tag = "SDX11DDVV08SCC"
    #cal_tag = "SDX11DDVVC"
    #cal_tag = "DX11DDVVSLOW"
    #cal_tag = "SWAT1C"
    #cal_tag = "DX11DDVVSLOW"
    #conf["input_map_filename"] = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx11_delta/lfi/LFI_SkyMap_%03d_1024_DX11D_full.fits" % freq

    # peter calibration
    conf["cal_tag"] = "CONVSOLDIP"
    conf["decalibrate_cal"] = None
    conf["input_map_filename"] = None
    conf["add_white_noise"] = False

    # open volts file
    efftype = "C"
    filename = latest_exchange(freq, ods=od, exchangefolder = "/project/projectdirs/planck/data/mission/lfi_ops_dx11_delta/", type = efftype)
    with pyfits.open(filename) as fitsfile:

        for ch in pl.f[freq].ch:
            # open poiting
            pnt = pointing.DiskPointing(od, ch.f.freq)
            ch_pnt = ch.inst[ch.tag.replace("S","M")]
            pix, psi = pnt.get_pix_psi(ch_pnt, NSIDE, nest=False)
            
            if ch.arm == "S":
                # S is using the same pointing as M, so also the same psi,
                # to recover the S psi, we remove the M PSI_POL angle (~90 deg) which
                # was added by the pointing code
                psi -= np.radians(ch_pnt.get_instrument_db_field("PSI_POL"))
                psi[psi < - np.pi] += 2*np.pi

            filetag = [ext.name for ext in fitsfile if ch.tag.replace("-","_").upper() in ext.name][0]
            obt = fitsfile["OBT"].data["OBT"]/2**16

            dipole = get_4pi_dipole(obt, pix, psi, ch, "solar_system")

            fitsfile[filetag].data[filetag][:] = dipole
            assert np.isnan(fitsfile[filetag].data[filetag]).sum() == 0

        formatted_date = '20110101'
        filename = '%s%03d-%04d-%s-%s.fits' % ("L", freq, od, "R", formatted_date)
        no_noise_tag = ""
        output_folder = "/global/scratch2/sd/planck/user/zonca/data/%s/%04d/" % (conf["cal_tag"] + no_noise_tag, od)
        try:
            os.mkdir(output_folder)
        except:
            pass
        fitsfile.writeto(output_folder + filename, clobber=True)

if __name__ == "__main__":
    import sys
    do1od(int(sys.argv[1]), int(sys.argv[2]))
