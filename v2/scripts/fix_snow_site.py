#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Cold-season snow top-up of one or more towers' forcing with the B-4-1 routine
(scripts/sites/fix_forcing_snow.py repair_snow, reused, not copied): for every UTC
month from October to April the shortfall of tower snow-equivalent precipitation
against WFDE5 v2.1 CRU snowfall of the nearest 0.5-degree cell is added at the
WFDE5 snowfall hours with tower Tair at or below 0 C; recorded precipitation is
never removed. Only the Met file's Precip of the named towers is rewritten (the
data set keeps one copy, data/FLUXNET-CH4/README.md section 6); the FI-Lom
rescaling of the original script is not touched. The monthly table goes to
v2/results/snow_topup_<sites>.csv.
Usage: fix_snow_site.py <tower,tower,...> [--apply]   (run on a compute node)"""
import sys

sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/sites')
import fix_forcing_snow as fs  # noqa: E402
import pandas as pd  # noqa: E402

V2 = '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2'


def main():
    sites = [s for s in sys.argv[1].split(',') if s]
    apply = '--apply' in sys.argv[2:]
    tabs = [fs.repair_snow(s, apply) for s in sites]
    out = f"{V2}/results/snow_topup_{'_'.join(sites)}.csv"
    pd.concat(tabs).to_csv(out, index=False)
    print(f'monthly table: {out}')
    if not apply:
        print('dry run: nothing written (add --apply)')


if __name__ == '__main__':
    main()
