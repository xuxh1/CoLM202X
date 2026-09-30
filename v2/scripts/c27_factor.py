"""Offline C-27 factor per tower (V-38 prediction): annual mean soil temperature
of the top metre from the run's history, and the production-weighted ratio of
the new temperature factor to the old one, 2**((295 - T_ann)/10).
Usage: c27_factor.py <case under cases/> <tag>"""
import glob, sys
import numpy as np, netCDF4 as nc
C = '/share/home/dq076/mode/Methane/cases/'
DZ = np.diff([0., 0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289])
case, tag = sys.argv[1:3]
for p in sorted(glob.glob(f'{C}{case}/sites_{tag}/*')):
    s = p.split('/')[-1]
    if s == '_conf':
        continue
    t = []
    for f in sorted(glob.glob(f'{p}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as d:
            a = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)[:, 0, 5:12]
            t.append((a * DZ).sum(1) / DZ.sum())
    if not t:
        continue
    t = np.concatenate(t)
    ta = np.nanmean(t)
    print(f'{s:7s} T_ann {ta - 273.15:5.1f} C  factor {2 ** ((295 - ta) / 10):4.2f}')
