#!/usr/bin/env python
import sys
import numpy as np
import os.path
import os
from glob import glob
import envoy

def get_freq_from_chtag(chtag):
    horn = int(chtag[0:3])
    if horn <= 23:
        return 70
    elif horn <= 26:
        return 44
    else:
        return 30

single_ch_scale_down = {30:4., 44:8., 70:8.}

tag = sys.argv[1].replace("runs_", "")

run_configuration = {
    30 : dict(
        full =   dict(minutes=12, cores=600/24/2),
        survey = dict(minutes=8, cores=120/24)
        ),
    44 : dict(
        full =   dict(minutes=12, cores=960/24/2),
        survey = dict(minutes=8, cores=480/24/2)
        ),
    70 : dict(
        full =   dict(minutes=14, cores=4800/24/2),
        survey = dict(minutes=8, cores=480/24/2)
        ),
}

known_map_types = "survey", "full", "year"

try:
    os.mkdir(os.environ["SCRATCH"] + "/maps/%s" % tag)
except:
    pass

with open("madam_template_py.pbs") as f:
    madam_template = f.read()

for run_fullpath in sorted(glob("runs_{}/dx11_*_escratch.bin".format(tag))):
    run=os.path.basename(run_fullpath)
    print(run)
    run_tag = run.split("_")[2]
    run_chtag = run.split("_")[3]
    print(run_tag)
    single_channel = run_chtag.endswith("M") or run_chtag.endswith("S")
    freq = int(run_tag)
        
    name=run.replace("_escratch.bin", "")
    job=name.replace("dx11_delta_", "")
    madam_map = glob(os.environ["SCRATCH"] + "/maps/%s/%s_map.fits" % (tag, name))

    # pbs configuration
    map_type = None
    for each in known_map_types:
        if each in name:
            map_type = each


    if map_type:
        
        # pbs config
        if map_type in ["full", "survey"]:
            minutes = run_configuration[freq][map_type]["minutes"]
            cores =   run_configuration[freq][map_type]["cores"]
        elif map_type in ["fullnos2s4"]:
            minutes = run_configuration[freq]["full"]["minutes"]
            cores =   run_configuration[freq]["full"]["cores"]
        elif map_type in ["year"]:
            if freq == 70:
                minutes = run_configuration[freq]["full"]["minutes"]/2
            else:
                minutes = run_configuration[freq]["full"]["minutes"]
            cores =   run_configuration[freq]["full"]["cores"]/2

        queue = "regular"
        if map_type == "survey" or freq == 30:
            queue = "regular"

        if single_channel:
            cores = int(np.ceil(cores / single_ch_scale_down[freq]))

        temperature_only = "T" if single_channel else "F"

        if madam_map:
            print(madam_map[0] + " already exists, SKIP " + name)
        else:
            print("SUBMIT " + name)
            pbs_command = madam_template.format(RUN=run, TAG=tag, NAME=name, JOB=job, FREQ=freq, MINUTES=minutes, CORES=int(cores), QUEUE=queue, temperature_only=temperature_only)
            cmd_script_path = run_fullpath.replace(".bin", ".cmd")
            with open(cmd_script_path, "w") as f:
                f.write(pbs_command)
            print(pbs_command)
            #sys.exit(0)
            r = envoy.run('sbatch ' + cmd_script_path, timeout=20)

    else:
            print("UNKNOWN MAP TYPE " + name)


print("Modified madam_lfi.par, set kfirst and kfilter to F")
