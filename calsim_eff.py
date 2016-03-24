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

def rescan_map(input_map_filename, pix, psi):
    input_map = hp.read_map(input_map_filename, (0,1,2))
    signal = input_map[0][pix]
    qw, uw = compute_pol_weigths(psi)
    signal += input_map[1][pix] * qw
    signal += input_map[2][pix] * uw
    return signal

class CalibrationProcessor:

    def __init__(self, od):
        conn = sqlite3.connect(private.database)
        self.c = conn.cursor()
        self.c.execute('select start_time, pointID_unique, end_time from list_ahf_infos where start_time >= 106743583006869 and end_time > start_time and od = %d'% od)
        od_rings = self.c.fetchall()
        self.pid_obts = {}
        self.pid_obts_end = {}
        for i, (start_time, pointID_unique, end_time) in enumerate(od_rings):
            pid = int(pointID_unique.split('-')[0])
            self.pid_obts[pid] = start_time/2**16
            self.pid_obts_end[pid] = end_time/2**16
         
        self.pid_obts = pd.DataFrame({"obt":self.pid_obts})

    def interpolate_calibration(self, cal, obt):
        obt_cal = cal.gain.reindex(self.pid_obts.index).fillna(method="ffill").fillna(method="bfill")
        if np.isnan(obt_cal).sum() > 0:
            obt_cal[:] = cal.gain.mean()
        return interpolate.interp1d(self.pid_obts.obt, obt_cal, kind="zero", bounds_error=False, fill_value=obt_cal.mean())(obt)

def do1od(od, freq):
#od = 100
#freq = 30

    pl = Planck()

    NSIDE = 1024

    def get_dipole(obt, pix, ch, orbtype):
        from dipole import Dipole
        dip = Dipole(obt=obt, type=orbtype) 
        return dip.get(ch, np.array(hp.pix2vec(NSIDE, pix, nest=False)).T)

    def get_4pi_dipole(obt, pix, psi, ch, orbtype):
        from dipole import Dipole
        dip = Dipole(obt=obt, type=orbtype) 
        theta, phi = hp.pix2ang(NSIDE, pix, nest=False)
        return dip.get_4piconv_dx10(ch, theta, phi, psi)

    print("-----OD %d" % od)

    calibration_processor = CalibrationProcessor(od)

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
    #conf["decalibrate_cal"] = "DX11D"
    #conf["cal_tag"]="SS_DX11D"
    #conf["cal_tag"]="S_DX11D_C_8"
    conf["straylight"] = True

    # peter calibration
    # conf["cal_tag"] = "DX11DDVVMODEL"
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

            dipole = get_4pi_dipole(obt, pix, psi, ch, "total")

            if conf["input_map_filename"]:
                fitsfile[filetag].data[filetag][:] = rescan_map(conf["input_map_filename"], pix, psi)
                assert np.isnan(fitsfile[filetag].data[filetag]).sum() == 0
            else:
                fitsfile[filetag].data[filetag][:] = 0

            if conf["straylight"]:
                straylight = pd.read_hdf("/global/project/projectdirs/planck/software/zonca/dev/chi2cal/ringsets/galactic_straylight_%s_%d.h5"  % (ch.tag, 128), "data")
                pix = pnt.get_pix(ch_pnt, 128, nest=True)
                for pid, pid_obt in calibration_processor.pid_obt.items():
                    pid_obt_end = calibration_processor.pid_obt_end[pid]
                    straylight_pid = straylight.loc[pid]
                    pix_pid = pix[pid_obt < obt < pid_obt_end]
                    straylight_timeline = pd.Series(np.array(straylight_pid.straylight), index=straylight_pid.pix).reindex(pix_pid).fillna(method="ffill").fillna(method="bfill")
                    fitsfile[filetag].data[filetag][pid_obt < obt < pid_obt_end] += straylight_timeline
                assert np.isnan(fitsfile[filetag].data[filetag]).sum() == 0
                    
            fitsfile[filetag].data[filetag] += dipole

            if conf["add_white_noise"]:
                white_noise_sigma = ch.white_noise_sigma
                fitsfile[filetag].data[filetag] += np.random.normal(size=len(obt)) * np.sqrt(white_noise_sigma)
                assert np.isnan(fitsfile[filetag].data[filetag]).sum() == 0

        #        #baselines_file ="/global/project/projectdirs/planck/data/mission/baselines/lfi/dx11_delta/base_dx11_delta_%03d_full.fits"
        #        #madam_baselines = MadamBaselines(baselines_file % ch.f.freq)
            #rings = pd.DataFrame(c.fetchall(), columns=["start_time", "pointID_unique", "end_time"])

            if conf["decalibrate_cal"]:
                cal = rings.load_fits_gains(conf["decalibrate_cal"], ch.tag, by_ring=True)
                fitsfile[filetag].data[filetag] /= calibration_processor.interpolate_calibration(cal, obt)

            # calibrate and dipole-remove
            if "cal_tag" in conf:
                cal = rings.load_fits_gains(conf["cal_tag"], ch.tag, by_ring=True)
                fitsfile[filetag].data[filetag] *= calibration_processor.interpolate_calibration(cal, obt)
                fitsfile[filetag].data[filetag] -= dipole
                assert np.isnan(fitsfile[filetag].data[filetag]).sum() == 0

        formatted_date = '20110101'
        filename = '%s%03d-%04d-%s-%s.fits' % ("L", freq, od, "R", formatted_date)
        no_noise_tag = "WN" if conf["add_white_noise"] else ""
        ctag = conf.get("cal_tag", "dipoles")
        output_folder = "/global/scratch2/sd/planck/user/zonca/data/%s/%04d/" % (ctag + no_noise_tag, od)
        try:
            os.mkdir(output_folder)
        except:
            pass
        fitsfile.writeto(output_folder + filename, clobber=True)

if __name__ == "__main__":
    import sys
    do1od(int(sys.argv[1]), int(sys.argv[2]))
