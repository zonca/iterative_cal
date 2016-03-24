from IPython.parallel import Client

def do1od(od, freq):
    import sys
    sys.path.append("/global/project/projectdirs/planck/software/zonca/dev/iterative_cal")
    import calsim_eff
    calsim_eff.do1od(od, freq)

rc = Client(profile="paral")

lview = rc.load_balanced_view()

od = list(range(933, 1540)) + list(range(1541, 1604+1))
od = list(range(153, 255))
od = list(range(91, 1540)) + list(range(1541, 1545+1))

#for freq in [30, 44,70]:
for freq in [70]:
    print(freq)
    lview.map(do1od, od, [freq]*len(od), block=False)
