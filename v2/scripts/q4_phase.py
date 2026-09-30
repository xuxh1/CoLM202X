#!/usr/bin/env python
"""Q-4: why the CH4 monthly correlation drops at some towers under C-11.
Monthly climatology of CH4 flux (mg m-2 d-1), inundated fraction and water
level (D-2) for V-1, V-2 and the outflow-only ablation, with the observed CH4 and water table."""
import glob
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, f'{R}/scripts/sites')
import score_series as SS                                           # noqa: E402
RUNS = (('V1', 'paper_v2/v260924a/sp_v2', 'v2_base_sp30'), ('V2', 'paper_v2/v260924b/sp_v2_c11', 'v2_c11_sp30'),
        ('V2a', 'paper_v2/v260924b/sp_v2_c11', 'v2_c11a_sp30'))   # V2a: outflow only, step F
OB = pd.read_csv(f'{R}/data/_tmp/分析_260923_r6/观测.csv')


def clim(case, tag, s):
    rows = []
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        y = int(f[-7:-3])
        with nc.Dataset(f) as h, nc.Dataset(f.replace('_hist_', '_hist_tracer_')) as e:
            g = lambda d, v: np.ma.filled(d[v][:], np.nan).ravel()
            c = g(e, 'f_methane_surf_flux_tot_active')
            if c.size != 12:
                continue
            fi, z, w = g(e, 'f_methane_finundated'), g(h, 'f_wetzwt'), g(h, 'f_wdsrf')
            cs, cu = g(e, 'f_methane_surf_flux_tot_sat'), g(e, 'f_methane_surf_flux_tot_unsat')
            for i in range(12):
                rows.append({'m': i + 1, 'ch4': c[i] * SS.MOD_TO_MGM2D, 'sat': cs[i] * SS.MOD_TO_MGM2D,
                             'unsat': cu[i] * SS.MOD_TO_MGM2D, 'F': fi[i],
                             'lev': -w[i] / 1000 if w[i] > 1 else z[i]})
    return pd.DataFrame(rows).groupby('m').mean()


args = sys.argv[1:]
if args and args[0] == '--val':                                    # validation tags of V-4
    RUNS = (('V1', 'paper_v2/v260924a/sp_v2', 'v2_base_val'), ('V2', 'paper_v2/v260924b/sp_v2_c11', 'v2_c11_val'))
    args = args[1:]
for s in args:
    obs = SS.obs_monthly(s)
    oc = pd.Series({k: v for k, v in obs.items()}).groupby(level=1).mean()
    ow = OB[OB.site == s].groupby('month').wtd_obs.mean()
    print(f'== {s}')
    t = pd.DataFrame({'obs_ch4': oc, 'obs_wt': ow})
    for lab, case, tag in RUNS:
        c = clim(case, tag, s)
        for k in ('ch4', 'sat', 'unsat', 'F', 'lev'):
            t[f'{lab}_{k}'] = c[k]
    print(t.round(2).T.to_string())
