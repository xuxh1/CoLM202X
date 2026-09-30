#!/usr/bin/env python
"""Step 2 premise: does the GLWD surface data's own soil-to-wetland area ratio
of each tower's 2-degree cell match the ratio the offline store needed
(round 3, inflow = ratio x upland subsurface drainage of the global run)?
Cell fractions from landdata/diag/patchfrac_elm_2005.nc of the V2 global
case (TypeIndex 11 wetland; 13 urban, 15 ice, 17 water excluded from soil)."""
import csv

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
PF = f'{R}/cases/paper_v2/v260924c/g2_v2/landdata/diag/patchfrac_elm_2005.nc'
SRF = f'{R}/cases/v260923c/sp_paper_wtdfloor/sites_r8_dyn_bf_sp30'
t3 = pd.read_csv(f'{R}/data/_tmp/分析_260924_连通离线/逐站_第三轮.csv').set_index('site')
with nc.Dataset(PF) as h:
    print({k: h[k].shape for k in h.variables})
    pf = np.ma.filled(h['patchfrac_elm'][:], 0.0)
    lat = h['lat'][:] if 'lat' in h.variables else None
    lon = h['lon'][:] if 'lon' in h.variables else None
rows = []
for s in t3.index:
    with nc.Dataset(f'{SRF}/{s}/landdata/srfdata.nc') as g:
        la, lo = float(g['latitude'][:]), float(g['longitude'][:])
    i = int(np.argmin(abs(lat - la))); j = int(np.argmin(abs(((lon - lo + 180) % 360) - 180)))
    f = pf[:, i, j]
    wet = f[11]; soil = f.sum() - f[11] - f[13] - f[15] - f[17] - f[0]
    rows.append({'site': s, 'cls': t3.loc[s, 'cls'], 'wet_frac': wet, 'soil_frac': soil,
                 'ratio_cell': soil / wet if wet > 0 else np.inf,
                 'ratio_needed_rsub': t3.loc[s, 'x_rsub'], 'ratio_needed_tau30': t3.loc[s, 'x_tau30']})
t = pd.DataFrame(rows)
print(t.round(3).to_string(index=False))
