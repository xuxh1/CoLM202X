#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Per-tower inputs of the finite lateral inflow (C-67) for a site config.
For every tower of <cfg>/LIST_sites.csv:
  DEF_WETLAND_INFLOW_SITE_RUNOFF(12)  upland runoff, mm/day, monthly
      climatology of the grid-mean total runoff f_rnof of the global O-5 run
      (cases/paper_v2/v260928b/g2_mr, before the floodplain re-infiltration C-50, whose river water cycled through the soil lifts the grid runoff of floodplain cells of later runs to 16-31 mm/day) in the
      2-degree cell holding the tower, floored at 0. The grid mean includes
      the cell's wetland tiles, whose floor refill books negative runoff, so
      in wetland-rich cells it is an underestimate of the soil tiles' runoff.
  DEF_WETLAND_INFLOW_RATIO_SITE        min(R, cap), R the upland-to-wetland
      area ratio of the HydroBASINS level-12 sub-basin holding the tower
      (v2/results/hybas_ratio.log); towers without an R take the cap, the
      value almost every 2-degree cell of the global run reaches.
The two columns are merged into <cfg>/SITE_MAIN.csv (rows added for towers
not yet in it); DEF_WETLAND_LATERAL_INFLOW itself is set in template.nml
(cfg_variant.py). Tower positions from the srfdata of a built site tag.
Usage: mk_inflow_site.py <cfg> <cap>   (run on a compute node)"""
import csv
import glob
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

V2 = Path('/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2')
G = '/share/home/dq076/mode/Methane/cases/paper_v2/v260928b/g2_mr/history'
SRF = '/share/home/dq076/mode/Methane/cases/paper_v2/v260929/sp_v2_c52/sites_v2_oo5'
K_RUN, K_R = 'DEF_WETLAND_INFLOW_SITE_RUNOFF', 'DEF_WETLAND_INFLOW_RATIO_SITE'


def hybas_ratio():
    out = {}
    for line in (V2 / 'results' / 'hybas_ratio.log').read_text().splitlines():
        p = line.split('|')
        if len(p) >= 5 and p[0].split() and '-' in p[0].split()[0] and p[3].split() and p[3].split()[0] != 'nan':
            try:
                out[p[0].split()[0]] = float(p[3].split()[0])
            except ValueError:
                pass
    return out


def runoff_clim():
    files = [f for f in sorted(glob.glob(f'{G}/g2_mr_hist_2*.nc'))
             if len(nc.Dataset(f)['time']) == 12]
    acc = None
    for f in files:
        with nc.Dataset(f) as h:
            v = np.ma.filled(h['f_rnof'][:].astype(float), np.nan)
            lat, lon = h['lat'][:], h['lon'][:]
        acc = v if acc is None else acc + v
    return acc / len(files), np.asarray(lat), np.asarray(lon), [Path(f).name for f in files]


def main():
    cfg, cap = V2 / 'config' / sys.argv[1], float(sys.argv[2])
    R = hybas_ratio()
    clim, lat, lon, used = runoff_clim()
    sites = [r['ID'] for r in csv.DictReader(open(cfg / 'LIST_sites.csv'))]
    vals = {}
    for s in sites:
        with nc.Dataset(f'{SRF}/{s}/landdata/srfdata.nc') as g:
            la, lo = float(g['latitude'][:]), float(g['longitude'][:])
        i = int(np.argmin(abs(lat - la)))
        j = int(np.argmin(abs(((lon - lo) + 180) % 360 - 180)))
        m = np.maximum(np.nan_to_num(clim[:, i, j], nan=0.0), 0.0) * 86400.
        r = min(R[s], cap) if s in R else cap
        vals[s] = (', '.join(f'{x:.3f}' for x in m), f'{r:.2f}')
        print(f'{s:8s} r {r:5.2f} (R {R.get(s, float("nan")):7.2f})  runoff mm/d '
              f'ann {m.mean():6.3f}  min {m.min():6.3f}  max {m.max():6.3f}')
    p = cfg / 'SITE_MAIN.csv'
    rows = list(csv.DictReader(open(p, newline=''))) if p.is_file() else []
    cols = list(rows[0].keys()) if rows else ['ID']
    for k in (K_RUN, K_R):
        if k not in cols:
            cols.append(k)
    have = {r['ID'] for r in rows}
    rows += [{'ID': s} for s in sites if s not in have]
    for r in rows:
        if r['ID'] in vals:
            r[K_RUN], r[K_R] = vals[r['ID']]
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, cols, restval='')
        w.writeheader()
        w.writerows(rows)
    print(f'runoff from {used}; wrote {len(vals)} towers into {p}')


if __name__ == '__main__':
    main()
