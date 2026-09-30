#!/usr/bin/env python
"""Q-3: month-by-month climatology of the wetland patch's assimilation
(f_assim, r9_floor_bf_sp30) against observed GPP, with LAI and 2 m air
temperature, for four towers from cold to warm."""
import glob

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
CASE = f'{R}/cases/v260924/sp_paper_frozenpsi/sites_r9_floor_bf_sp30'
OBS = f'{R}/data/FLUXNET-CH4/Observation'
for s in ('SE-Deg', 'FI-Sii', 'NZ-Kop', 'US-Tw1'):
    rows = []
    for f in sorted(glob.glob(f'{CASE}/{s}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as h:
            g = lambda v: np.ma.filled(h[v][:], np.nan).ravel()
            if g('f_assim').size != 12:
                continue
            for i in range(12):
                rows.append({'m': i + 1, 'assim': g('f_assim')[i] * 1e6, 'lai': g('f_lai')[i], 'tref': g('f_tref')[i] - 273.15})
    m = pd.DataFrame(rows).groupby('m').mean()
    fo = sorted(glob.glob(f'{OBS}/{s}_*_Flux.nc'))[0]
    with nc.Dataset(fo) as d:
        t = nc.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        o = pd.Series(np.ma.filled(d['GPP'][:], np.nan).ravel().astype(float), index=pd.DatetimeIndex(t))
    m['gpp_obs'] = o.groupby(o.index.month).mean()
    print(s); print(m.round(2).T.to_string())
