#!/usr/bin/env python
"""Is the wetland patch's biophysical photosynthesis (f_assim, prescribed LAI)
a usable plant-carbon signal? Compare with observed GPP per tower over matched
months (monthly means need >= 480 valid half-hours): annual totals, ratio,
cross-tower rank correlation, within-tower monthly correlation.
Usage: gpp_check.py [<case under cases/> <tag>]  (default: the r9_floor_bf_sp30 tag)"""
import csv
import glob
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
CASE = (f'{R}/cases/{sys.argv[1]}/sites_{sys.argv[2]}' if len(sys.argv) > 2
        else f'{R}/cases/v260924/sp_paper_frozenpsi/sites_r9_floor_bf_sp30')
OBS = f'{R}/data/FLUXNET-CH4/Observation'
SITES = [r['ID'] for r in csv.DictReader(open(f'{R}/scripts/sites/LIST_sites_wtd23.csv'))]
SALT = {'DE-Hte', 'US-Srr', 'US-StJ'}
rows = []
for s in SITES:
    if s in SALT:
        continue
    m = {}
    for f in sorted(glob.glob(f'{CASE}/{s}/history/{s}_hist_2*.nc')):
        y = int(f[-7:-3])
        with nc.Dataset(f) as h:
            a = np.ma.filled(h['f_assim'][:], np.nan).ravel()
            if a.size != 12:
                continue
            for i in range(12):
                m[pd.Timestamp(y, i + 1, 1)] = a[i]
    m = pd.Series(m, name='assim') * 1e6              # umol m-2 s-1
    f = sorted(glob.glob(f'{OBS}/{s}_*_Flux.nc'))[0]
    with nc.Dataset(f) as d:
        t = nc.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        g = np.ma.filled(d['GPP'][:], np.nan).ravel().astype(float) if 'GPP' in d.variables else np.full(len(t), np.nan)
    o = pd.Series(g, index=pd.DatetimeIndex(t))
    cnt, mean = o.resample('MS').count(), o.resample('MS').mean()
    o = mean.where(cnt >= 480).rename('gpp')
    j = pd.concat([m, o], axis=1).dropna()
    if len(j) < 6:
        rows.append({'site': s, 'n': len(j)})
        continue
    to_g = 12.011e-6 * 86400 * 365.25
    rows.append({'site': s, 'n': len(j), 'assim_ann': j.assim.mean() * to_g, 'gpp_obs_ann': j.gpp.mean() * to_g,
                 'ratio': j.assim.mean() / j.gpp.mean(), 'r_month': j.assim.corr(j.gpp)})
t = pd.DataFrame(rows)
print(t.round(2).to_string(index=False))
x = t.dropna(subset=['ratio'])
print('towers %d; median ratio model/obs %.2f; spearman of annual %.2f; median monthly r %.2f' % (
    len(x), x.ratio.median(), x.assim_ann.corr(x.gpp_obs_ann, method='spearman'), x.r_month.median()))
