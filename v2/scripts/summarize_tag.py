#!/usr/bin/env python
"""Score one V2 site tag and write the summary into v2/results/.
1. runs v2/scripts/score_series_v2.py (score_series with the D-17 gate, 23-tower list) and moves its txt/csv
   from the main tree's outputs/analysis into v2/results/scores/;
2. CH4: KGE_ln, beta, alpha, r medians over the 20 non-saline towers, per
   class, and for all 23;
3. water level (ledger D-2): the model level is minus the ponded depth when
   the pond exceeds 1 mm, else the water-table depth; per tower the median
   over months with an observed depth and the monthly correlation.
Usage: summarize_tag.py <case under cases/> <tag>"""
import csv
import glob
import os
import shutil
import subprocess
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
V2 = f'{R}/CoLM202X-paper-v2/v2'
PY = '/share/home/dq076/software/miniconda3/envs/py311/bin/python'
SITES = {r['ID']: int(r['SITE_wetland_class']) for r in csv.DictReader(open(f'{R}/scripts/sites/LIST_sites_wtd23.csv'))}
CLS = {1: 'permafrost', 3: 'bog', 4: 'fen', 5: 'marsh', 6: 'salt_marsh', 7: 'trop_swamp'}
SALT = {'DE-Hte', 'US-Srr', 'US-StJ'}
OB = pd.read_csv(f'{R}/data/_tmp/分析_260923_r6/观测.csv')


def water_level(case, tag, site):
    rows = []
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{site}/history/{site}_hist_2*.nc')):
        y = int(f[-7:-3])
        with nc.Dataset(f) as h:
            z = np.ma.filled(h['f_wetzwt'][:], np.nan).ravel()
            w = np.ma.filled(h['f_wdsrf'][:], np.nan).ravel()
        if z.size != 12:
            continue
        for i in range(12):
            rows.append({'year': y, 'month': i + 1, 'lev': -w[i] / 1000.0 if w[i] > 1.0 else z[i]})
    return pd.DataFrame(rows)


def main():
    case, tag = sys.argv[1], sys.argv[2]
    subprocess.run([PY, '-W', 'ignore', f'{V2}/scripts/score_series_v2.py', tag, f'{R}/cases/{case}',
                    '--sitelist', f'{R}/scripts/sites/LIST_sites_wtd23.csv'], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.makedirs(f'{V2}/results/scores', exist_ok=True)
    for ext in ('txt', 'csv'):
        shutil.move(f'{R}/outputs/analysis/score_series_{tag}.{ext}', f'{V2}/results/scores/score_series_{tag}.{ext}')
    s = pd.read_csv(f'{V2}/results/scores/score_series_{tag}.csv')
    s['cls'] = s.site.map(lambda x: CLS[SITES[x]])
    L = [f'# {tag} ({case}); CH4 medians KGE_ln / beta / alpha / r']
    four = lambda d: ' / '.join(f'{d[k].median():.2f}' for k in ('KGEln', 'beta', 'alpha', 'r'))
    L.append(f'20 non-saline: {four(s[~s.site.isin(SALT)])}')
    L.append(f'all 23:        {four(s)}')
    for c in ('permafrost', 'bog', 'fen', 'marsh', 'salt_marsh', 'trop_swamp'):
        d = s[s.cls == c]
        L.append(f'  {c:11s} n={len(d):2d}  {four(d)}')
    wl = []
    for site in SITES:
        if site in SALT:
            continue
        m = water_level(case, tag, site)
        j = m.merge(OB[OB.site == site][['year', 'month', 'wtd_obs']], on=['year', 'month']).dropna()
        if len(j) < 6:
            continue
        r = j.lev.corr(j.wtd_obs) if j.lev.std() > 1e-6 else np.nan
        wl.append({'site': site, 'obs': j.wtd_obs.median(), 'mod': j.lev.median(),
                   'err': abs(j.lev.median() - j.wtd_obs.median()), 'r': r})
    w = pd.DataFrame(wl)
    L.append(f'water level (D-2), {len(w)} towers: median |median error| {w.err.median():.2f} m, median r {w.r.median():.2f}')
    L.append(w.round(2).to_string(index=False))
    L.append('')
    L.append(s[['site', 'cls', 'KGEln', 'beta', 'alpha', 'r']].round(2).to_string(index=False))
    out = f'{V2}/results/{tag}_summary.txt'
    open(out, 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:10]))


if __name__ == '__main__':
    main()
