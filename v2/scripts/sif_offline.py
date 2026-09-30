#!/usr/bin/env python
"""Offline estimate for C-14: the CLM4Me seasonal inundation factor
S = (beta (f - fbar) + fbar) / f, S <= 1 (Riley et al. 2011, appendix B;
CLM5 tech note eq. 24.21), applied to the saturated part of the soil tile
(floodplain). f = monthly soil-tile inundated fraction, fbar = its annual mean
weighted by heterotrophic respiration (grid f_hr as the weight), beta = 0.2.
Emission of each floodplain cell-month is scaled by S (production in the
saturated column scales by S; transport and oxidation assumed proportional).
V-3 years 2000-2001, by latitude band, with per-area values."""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
BETA = 0.2
c = B.Case('paper_v2/v260924c', 'g2_v2')
A_land, _, _ = VB.load_areas(c, 2000)
BANDS = (('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90), ('90S-30S', -90, -30), ('all', -90, 90))
acc = {}
for y in (2000, 2001):
    tr = xr.open_dataset(f'{c.dir}/history/g2_v2_hist_tracer_{y}.nc')
    mn = xr.open_dataset(f'{c.dir}/history/g2_v2_hist_{y}.nc')
    w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
    f = tr['f_methane_soil_finundated'].fillna(0)
    hr = mn['f_hr'].fillna(0).clip(min=0)
    fbar = (f * hr).sum('time') / hr.sum('time').where(hr.sum('time') > 0)
    S = ((BETA * (f - fbar) + fbar) / f.where(f > 0)).clip(max=1).fillna(1)
    e = tr['f_methane_surf_flux_soil'].fillna(0) * A_land * w * B.M_CH4
    As = tr['f_methane_area_soil'].fillna(0) * A_land
    d = {'E': e.sum('time'), 'E_S': (e * S).sum('time'), 'Amean': (f * As).mean('time'), 'Amax': (f * As).max('time'),
         'Sw': (S * e).sum('time')}
    for k, v in d.items():
        acc[k] = acc.get(k, 0) + v.load() / 2
    tr.close(); mn.close()
print('# floodplain (soil tile), 2000-2001; E Tg/yr; per-area g CH4 m-2 yr-1; S = emission-weighted mean factor')
for b, lo, hi in BANDS:
    m = (acc['E'].lat >= lo) & (acc['E'].lat < hi) & (acc['Amax'] > 0)
    E, ES = float(acc['E'].where(m).sum()), float(acc['E_S'].where(m).sum())
    Am, Ax = float(acc['Amean'].where(m).sum()) / 1e12, float(acc['Amax'].where(m).sum()) / 1e12
    print(f'{b:8s} E {E:6.1f} -> {ES:6.1f} (x{ES / E if E else np.nan:4.2f})  per mean flooded area {E / Am if Am else np.nan:6.1f} -> '
          f'{ES / Am if Am else np.nan:6.1f}  per max flooded area {E / Ax if Ax else np.nan:6.1f} -> {ES / Ax if Ax else np.nan:6.1f}')
