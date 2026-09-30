#!/usr/bin/env python
"""Q-1: which WAD2M v2.0 statistic gives the 7.9 (annual average) and 13.6
(peak season maximum) Mkm2 quoted by Saunois et al. (2025)? Candidates:
monthly-mean area; per-year sum of per-pixel annual maxima; per-pixel maximum
over the whole record; per-year maximum of the monthly global total."""
import netCDF4 as nc
import numpy as np

F = '/share/home/dq076/data/methane/WAD2M/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc'
with nc.Dataset(F) as h:
    print({k: h[k].shape for k in h.variables})
    lat = h['lat'][:] if 'lat' in h.variables else h['latitude'][:]
    fw = h['Fw']
    nt = fw.shape[0]
    R = 6371.0
    dlat = abs(float(lat[1] - lat[0]))
    area = (np.radians(dlat) ** 2) * R ** 2 * np.cos(np.radians(np.asarray(lat)))   # km2 per 0.25 cell (square in lon/lat)
    area2d = area[:, None]
    ann_max_sum, monthly_tot, overall_max = [], [], None
    for y in range(nt // 12):
        x = np.ma.filled(fw[y * 12:(y + 1) * 12].astype(np.float64), 0.0)
        x = np.clip(x, 0, 1)
        monthly_tot += list((x * area2d).sum(axis=(1, 2)) / 1e6)
        pm = x.max(axis=0)
        ann_max_sum.append((pm * area2d).sum() / 1e6)
        overall_max = pm if overall_max is None else np.maximum(overall_max, pm)
    monthly_tot = np.array(monthly_tot)
    print('years', nt // 12)
    print('mean of monthly totals          %.2f Mkm2' % monthly_tot.mean())
    print('mean of annual max of monthly   %.2f' % monthly_tot.reshape(-1, 12).max(1).mean())
    print('mean of per-pixel annual max    %.2f (range %.2f-%.2f)' % (np.mean(ann_max_sum), min(ann_max_sum), max(ann_max_sum)))
    print('per-pixel max over all years    %.2f' % ((overall_max * area2d).sum() / 1e6))
