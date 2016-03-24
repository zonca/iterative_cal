
from IPython.parallel import Client
rc = Client(profile="paral")
lview = rc.load_balanced_view()
def do1od(od, freq):
        import sys
        sys.path.append("/global/project/projectdirs/planck/software/zonca/dev/iterative_cal")
        from sol_dip_only_calsim_eff import do1od
        do1od(od, freq)
    
od = list(range(91, 1540)) + list(range(1541, 1544+1))
for freq in [30, 44,70]:
        print(freq)
        lview.map(do1od, od, [freq]*len(od), block=False)
    
