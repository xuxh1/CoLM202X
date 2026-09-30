#!/usr/bin/env python
"""Q-4 / C-15: monthly climatology at cold peatland towers for two tags:
CH4 flux (mg m-2 d-1), inundated fraction, water level (D-2), soil temperature
of the top three layers (C), CH4 production (mg m-2 d-1), with observed CH4.
Usage: spring_pulse.py <site> [<site> ...]"""
import glob
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, f'{R}/scripts/sites')
import score_series as SS                                           # noqa: E402
RUNS = (('C12a', 'paper_v2/v260924d/sp_v2_c12', 'v2_c12a_sp30'), ('C15', 'paper_v2/v260925d/sp_v2_c15', 'v2_c15_sp30'))


def clim(case, tag, s):
    rows = []
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as h, nc.Dataset(f.replace('_hist_', '_hist_tracer_')) as e:
            g = lambda d, v: np.ma.filled(d[v][:], np.nan)
            c = g(e, 'f_methane_surf_flux_tot_active').ravel()
            if c.size != 12:
                continue
            T = g(h, 'f_t_soisno').reshape(12, -1)
            ts = T[:, -10:][:, :3].mean(axis=1) - 273.15
            p = g(e, 'f_methane_prod_tot').ravel()
            fi, z, w = g(e, 'f_methane_finundated').ravel(), g(h, 'f_wetzwt').ravel(), g(h, 'f_wdsrf').ravel()
            for i in range(12):
                rows.append({'m': i + 1, 'ch4': c[i] * SS.MOD_TO_MGM2D, 'prod': p[i] * SS.MOD_TO_MGM2D, 'F': fi[i],
                             'lev': -w[i] / 1000 if w[i] > 1 else z[i], 'Ttop': ts[i]})
    return pd.DataFrame(rows).groupby('m').mean()


for s in sys.argv[1:]:
    obs = SS.obs_monthly(s)
    t = pd.DataFrame({'obs_ch4': pd.Series(obs).groupby(level=1).mean()})
    for lab, case, tag in RUNS:
        c = clim(case, tag, s)
        for k in ('ch4', 'prod', 'F', 'lev', 'Ttop'):
            t[f'{lab}_{k}'] = c[k]
    print(f'== {s}')
    print(t.round(2).T.to_string())
