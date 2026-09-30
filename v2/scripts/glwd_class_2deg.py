#!/usr/bin/env python
"""GLWD v2 wetland-tile classes (16-19, 22-27) aggregated to the 2-degree model
grid: area of each class per cell (km2) from the per-class percent layers
(pct/100 x pixel area), and the forested share of the cell's tile classes
(16, 18, 22, 24, 26 forested; 17, 19, 23, 25, 27 non-forested). Input for
C-13 (vegetation of the wetland tile by forested / non-forested share).
Output: v2/data/glwd_class_2deg.nc, rows from 84N (GLWD coverage 84N-56S)."""
import os

import netCDF4 as nc
import numpy as np
import rasterio

B = '/share/home/dq076/mode/Methane/data/GLWD_v2/area_by_class_pct/GLWD_v2_0_area_by_class_pct/'
OUT = '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/data/glwd_class_2deg.nc'
SET = list(range(16, 20)) + list(range(22, 28))
FOREST = [16, 18, 22, 24, 26]
R = 6371.0072
d = 1.0 / 240
NROW = 33600
step = 480                                           # one 2-degree row
src = {k: rasterio.open(B + f'GLWD_v2_0_class_{k:02d}_pct.tif') for k in SET}
cell = {k: np.zeros((NROW // step, 180)) for k in SET}
for i0 in range(0, NROW, step):
    w = rasterio.windows.Window(0, i0, 86400, step)
    lat_n = 84 - i0 * d - np.arange(step) * d
    area = (R ** 2 * np.radians(d) * (np.sin(np.radians(lat_n)) - np.sin(np.radians(lat_n - d))))[:, None]
    for k in SET:
        p = src[k].read(1, window=w).astype(np.float64)
        p[p > 100] = 0
        a = area * p / 100.0
        cell[k][i0 // step, :] = a.reshape(step, 180, 480).sum(axis=(0, 2))
    print('row', i0 // step, flush=True)
tot = sum(cell[k] for k in SET)
forest = sum(cell[k] for k in FOREST)
o = nc.Dataset(OUT, 'w')
o.createDimension('lat', tot.shape[0]); o.createDimension('lon', 180)
o.createVariable('lat', 'f8', ('lat',))[:] = 84 - 1 - 2 * np.arange(tot.shape[0])
o.createVariable('lon', 'f8', ('lon',))[:] = -180 + 1 + 2 * np.arange(180)
for k in SET:
    v = o.createVariable(f'area_class_{k:02d}', 'f8', ('lat', 'lon')); v[:] = cell[k]; v.units = 'km2'
v = o.createVariable('area_tile_classes', 'f8', ('lat', 'lon')); v[:] = tot; v.units = 'km2'
v = o.createVariable('forested_share', 'f8', ('lat', 'lon'))
v[:] = np.where(tot > 0, forest / np.where(tot > 0, tot, 1), np.nan)
v.long_name = 'forested classes (16,18,22,24,26) / all tile classes (16-19,22-27)'
o.source = 'GLWD v2.0 per-class percent layers; script v2/scripts/glwd_class_2deg.py'
o.close()
print('total tile-class area %.3f Mkm2, forested %.3f' % (tot.sum() / 1e6, forest.sum() / 1e6))
