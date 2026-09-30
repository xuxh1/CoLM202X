#!/usr/bin/env python
"""WAD2M v2.0 wetland area by latitude band (2010-2019): mean of per-pixel
annual maxima and mean of monthly area; the per-area emission the GCP band
targets imply (Zhang et al. 2025 table 1 ensemble means, v2/目标与判据.md);
and the model's inundated area by band (wetland tile + floodplain, per-cell
annual max and monthly mean) for V-3, years 2000-2001."""
import os
import sys

import netCDF4 as nc
import numpy as np
import xarray as xr

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B

F = '/share/home/dq076/data/methane/WAD2M/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc'
BANDS = (('90S-30S', -90, -30, 3.0), ('30S-30N', -30, 30, 110.0), ('30N-60N', 30, 60, 34.5), ('60N-90N', 60, 90, 12.0))
with nc.Dataset(F) as h:
    lat = np.asarray(h['lat'][:])
    dl = abs(float(lat[1] - lat[0]))
    area = (np.radians(dl) ** 2) * 6371.0 ** 2 * np.cos(np.radians(lat))      # km2 per cell
    amax = {b: [] for b, *_ in BANDS}
    amean = {b: [] for b, *_ in BANDS}
    for y in range(10, 20):                                                    # 2010-2019
        x = np.clip(np.ma.filled(h['Fw'][y * 12:(y + 1) * 12].astype(np.float64), 0.0), 0, 1)
        pm, mm = x.max(axis=0), x.mean(axis=0)
        for b, lo, hi, _ in BANDS:
            m = (lat >= lo) & (lat < hi)
            amax[b].append((pm[m] * area[m, None]).sum() / 1e6)
            amean[b].append((mm[m] * area[m, None]).sum() / 1e6)
c = B.Case('paper_v2/v260924c', 'g2_v2')
_, A_act, _ = VB.load_areas(c, 2000)
mod = {}
for yy in (2000, 2001):
    tr = xr.open_dataset(f'{c.dir}/history/g2_v2_hist_tracer_{yy}.nc')
    fa = tr['f_methane_finundated'].fillna(0) * A_act
    for b, lo, hi, _ in BANDS:
        s = fa.where((fa.lat >= lo) & (fa.lat < hi))
        mod.setdefault(b, []).append((float(s.max('time').sum()) / 1e12, float(s.mean('time').sum()) / 1e12))
    tr.close()
print('# band | WAD2M per-pixel annual max, monthly mean (Mkm2, 2010-2019) | GCP band target Tg/yr | implied g CH4 m-2 yr-1 per max area, per mean area | V-3 area annual max, monthly mean')
for b, lo, hi, e in BANDS:
    x, mn = np.mean(amax[b]), np.mean(amean[b])
    mx = np.mean([v[0] for v in mod[b]]), np.mean([v[1] for v in mod[b]])
    print(f'{b:8s} WAD2M {x:5.2f} {mn:5.2f} | target {e:6.1f} | implied {e / x:6.1f} {e / mn:6.1f} | V-3 {mx[0]:5.2f} {mx[1]:5.2f}')
