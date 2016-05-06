from planck import toast
import os
import sys
sys.path.append("/global/project/projectdirs/planck/software/keskital/python/general")

tag = sys.argv[1]

import planck_util as pu
import numpy as np

individual_detectors = False

#db = '/project/projectdirs/planck/data/mission/SIAM/rings_dx10.db'
#db = '/project/projectdirs/planck/data/mission/SIAM/rings_dx11.db'
db = '/global/homes/k/keskital/dx11d/aux/rings_dx11_delta.db'

#private.bad_ods = np.append(private.bad_ods, np.arange(454, 456)) # SCS switchover
#private.bad_ods = np.append(private.bad_ods, np.arange(938, 948)) # spin-up campaign
#private.bad_ods = np.append(private.bad_ods, np.arange(984, 985)) # Loss of thermal control, LFI reboot
#private.bad_rings = np.append(private.bad_rings, [
#        '00004200', # HFI data has NaNs on the first ring ...
#        '01050950', '01050960', '01050990', '01051000', # Unexpected STR switchover # 2 OD 242
#        '02066250', '02066410', # Unexpected STR switchover # 3 OD 288, 289
#        '93000720', # STR switchover test OD 560
#    ])

tag = sys.argv[1]

try:
    os.mkdir("runs_{}".format(tag))
except:
    pass


#for setname in ['030']:
for setname in ['030', '044', '070']:
    chtag = None

#for freqtag in ['030', '044', '070']:
#  for chtag in pu.list_planck(freqtag):
#    setname = freqtag + "_" + chtag.replace("LFI", "")

    subsets = {
        '030' : 30,
        '044' : 44,
        '070' : 70,
        #'070_18_23' : ['LFI18M', 'LFI18S', 'LFI23M', 'LFI23S'],
        #'070_19_22' : ['LFI19M', 'LFI19S', 'LFI22M', 'LFI22S'],
        #'070_20_21' : ['LFI20M', 'LFI20S', 'LFI21M', 'LFI21S'],
        '100' : 100,
        '100_1_4' : ['100-1a', '100-1b', '100-4a', '100-4b'],
        '100_2_3' : ['100-2a', '100-2b', '100-3a', '100-3b'],
        '143' : 143,
        '143_1_3' : ['143-1a', '143-1b', '143-3a', '143-3b'],
        '143_2_4' : ['143-2a', '143-2b', '143-4a', '143-4b'],
        '143_5' : ['143-5'],
        '143_6' : ['143-6'],
        '143_7' : ['143-7'],
        '217' : 217,
        '217_1' : ['217-1'],
        '217_2' : ['217-2'],
        '217_3' : ['217-3'],
        '217_4' : ['217-4'],
        '217_5_7' : ['217-5a', '217-5b', '217-7a', '217-7b'],
        '217_6_8' : ['217-6a', '217-6b', '217-8a', '217-8b'],
        '353' : 353,
        '353_1' : ['353-1'],
        '353_2' : ['353-2'],
        '353_3_5' : ['353-3a', '353-3b', '353-5a', '353-5b'],
        '353_4_6' : ['353-4a', '353-4b', '353-6a', '353-6b'],
        '353_7' : ['353-7'],
        '353_8' : ['353-8'],
        '545' : 545,
        '545_1' : ['545-1'],
        '545_2' : ['545-2'],
        '545_4' : ['545-4'],
        '857' : 857,
        '857_1' : ['857-1'],
        '857_2' : ['857-2'],
        '857_3' : ['857-3'],
        '857_4' : ['857-4'],
        }
    if chtag:
        subsets[setname] = [chtag]

    surveys = {
        'survey1' : [    3,  5483],
        'survey2' : [ 5484, 10957],
        'survey3' : [10958, 16454],
        'survey4' : [16455, 21482],
        'survey5' : [21483, 27404],
        'survey6' : [27405, 32761],
        'survey7' : [32762, 38574],
        'survey8' : [38575, 44072],
       # 'survey8' : [38575, 43711], # cut 1533 and after
       # #'survey9' : [44073, 45942],
        'year1'   : [    3, 10957],
        'year2'   : [10958, 21482],
        'year3'   : [21483, 32761],
        'year4'   : [32762, 43711], # cut 1533 and after
       # #'year13'  : [[    3, 10957],[21483, 32761]],
       # #'year24'  : [[10958, 21482],[32762, 44072]],
       # #'fullnosurvey24' : [[    3,  5483],[10958, 16454],[21483, 44072]],
       # #'year12'  : [    3, 21482],
       #'full'    : [    3, 44072], # survey9 is not included in LFI delivery
        'full'    : [    3, 43711], # cut 1533 and after
        'fullnos2s4'    : [    3, 43711], # cut 1533 and after
        #'hm1'     : [    3, 13234], # HFI half mission 1
        #'hm2'     : [13235, 27404], # HFI half mission 2
        #'full+'    : [    3, 45942], # This definition includes survey 9
        #'Bruce'   : [40552, 42891],
        }

    for subsetname, subset in subsets.items():

        if setname and setname not in subsetname:
            continue

        deaberrate = True
        extend_857 = True
        ahf_folder = '/global/cscratch1/sd/planck/data/mission/AHF_v7/'
        noise_tod = False
        eff_is_for_flags = False
        include_repointings = False

        ptcorfile = '/project/projectdirs/planck/data/mission/SIAM/PTCOR_REBASUN_1min.csv'
        ptcorfile = False
        wobble_high = False
        no_wobble = True

        #if int(subsetname[0:3]) < 100:
        fpdb = '/project/projectdirs/planck/data/mission/SIAM/LFI_RIMO_DX11_R2.0_M_arm_pnt.fits'
        calibration_file = None
        flagmask = 255
        dipole_removal = False
        flag_HFI_bad_rings = False
        nside = 1024
        psd = None
        obtmask = 1+4
        exchange_folder = '/project/projectdirs/planck/data/mission/lfi_ops_dx11_delta'
        exchange_folder = "/global/scratch2/sd/planck/user/zonca/data/"

        exchange_folder = os.environ["SCRATCH"] + "/data/"
        #exchange_folder += "S_DX11D_C_3"
        exchange_folder += tag
        efftype = 'R'
        #else:
        #    sfreq = setname[:3]
        #    ifreq = int( sfreq )
        #    fpdb = '/project/projectdirs/planck/data/mission/SIAM/HFI-RIMO-4_13_dx11.fits'
        #    #fpdb = '/global/homes/k/keskital/dx11d/aux/HFI-RIMO-4_13_dx11_v3.fits'
        #    calibration_file = None
        #    obtmask = 1 + 2 # unstable pointing + HCM flags
        #    flagmask = 1 # Just check the TotalFlag
        #    dipole_removal = False
        #    nside = 2048
        #    flag_HFI_bad_rings = '/project/projectdirs/planck/data/mission/SIAM/hfi_discarded_rings_dx11.txt'
        #    psd = None # '/project/projectdirs/planck/data/ffp6/mission/noise/noisefit_CHANNEL_ffp6.psd'
        #    efftype = 'R'
        #    exchange_folder = '/project/projectdirs/planck/data/mission/hfi_dx11'

        for survey, survey_span in surveys.items():

            if int(subsetname[0:3]) > 70: # HFI mission ends on LFI ring 27005-237
                if '13' in survey or '24' in survey or 'year3' in survey or 'fullnosurvey24' in survey: continue
                if survey_span[0] > 27005-237: continue
                if survey_span[1] > 27005-237: survey_span[1] = 27005-237
            else:
                if 'hm' in survey: continue # HFI definition of half mission

            if individual_detectors:
                detectors = pu.list_planck( subset )
                names = detectors
            else:
                detectors = [subset]
                names = [subsetname]

            for name, detector in zip(names, detectors):

                for platform in ['escratch',]:
                #for platform in ['project','hscratch','escratch']:
                    if int(subsetname[0:3]) > 70:
                        xmlfile = 'runs_' + tag + '/dx11_' + name + '_' + survey + '_' + platform + '.xml'
                    else:
                        xmlfile = 'runs_' + tag + '/dx11_delta_' + name + '_' + survey + '_' + platform + '.xml'

                    if os.path.isfile(xmlfile):
                        print 'Skipping existing file: ' + xmlfile
                        continue

                    if platform == 'hscratch':
                        remote_ahf_folder = ahf_folder.replace('project', 'scratch')
                        remote_exchange_folder = exchange_folder.replace('project', 'scratch')
                    elif platform == 'escratch':
                        remote_ahf_folder = ahf_folder.replace('project/projectdirs', 'scratch3/scratchdirs')
                        #remote_exchange_folder = exchange_folder.replace('project/projectdirs', 'scratch3/scratchdirs')
                        remote_exchange_folder = None                    
                    else:
                        remote_ahf_folder = None
                        remote_exchange_folder = None                    

                    toast_config = toast.ToastConfig(
                        db, fpdb, 
                        lfi_ring_range=survey_span, channels=detector, nside=nside,
                        ordering='NEST', coord='G', 
                        use_OCM=False,
                        calibration_file=calibration_file,
                        obtmask=obtmask, 
                        flagmask=flagmask, 
                        dipole_removal=dipole_removal,
                        exchange_folder=exchange_folder, 
                        remote_exchange_folder=remote_exchange_folder,
                        output_xml=xmlfile, 
                        ahf_folder=ahf_folder,
                        remote_ahf_folder=remote_ahf_folder,
                        flag_HFI_bad_rings=flag_HFI_bad_rings,
                        use_HCM=include_repointings,
                        ptcorfile=ptcorfile,
                        wobble_high=wobble_high,
                        no_wobble=no_wobble,
                        deaberrate=deaberrate,
                        extend_857=extend_857,
                        psd=psd,
                        noise_tod=noise_tod,
                        eff_is_for_flags=eff_is_for_flags,
                        efftype=efftype,
                        )

                    if int(subsetname[0:3]) < 100:
                        toast_config.ringdb.exclude_od( [ 455, 456, 1540 ] )

                    if survey.endswith("nos2s4"):
                        toast_config.ringdb.exclude_od( np.arange(270, 455) )
                        toast_config.ringdb.exclude_od( np.arange(606, 807) )
                        
                    # Exclude the spin-up campaign.

                    toast_config.ringdb.exclude_od( np.arange(939, 946) )
                    toast_config.ringdb.exclude_od( np.arange(1008, 1008+1) )
                    toast_config.run()
