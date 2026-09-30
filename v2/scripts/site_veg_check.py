"""Vegetation state of one site in a tag (main-set diagnosis): site-file land
cover and LAI, and the run's annual means of LAI, assimilation, PFT weights."""
import glob, sys
import numpy as np, netCDF4 as nc
C = '/share/home/dq076/mode/Methane/cases/'
case, tag, sid = sys.argv[1:4]
p = glob.glob(f'/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Sitedata/{sid}_*_site.nc')[0]
with nc.Dataset(p) as d:
    print(p.split('/')[-1])
    for v in d.variables:
        x = d[v][:]
        if np.size(x) <= 20:
            print(f'  {v}: {np.ravel(np.ma.filled(x, np.nan))}')
        else:
            a = np.ma.filled(x.astype(float), np.nan) if x.dtype.kind in 'fi' else None
            print(f'  {v}: shape {x.shape}' + (f' mean {np.nanmean(a):.3g} max {np.nanmax(a):.3g}' if a is not None else ''))
f = sorted(glob.glob(f'{C}{case}/sites_{tag}/{sid}/history/{sid}_hist_2*.nc'))[0]
with nc.Dataset(f) as d:
    print(f.split('/')[-1])
    for v in ('f_lai', 'f_sai', 'f_assim', 'f_fevpg', 'f_etr', 'f_zwt', 'f_wdsrf', 'f_fsno'):
        if v in d.variables:
            print(f'  {v}: {np.round(np.ma.filled(d[v][:].astype(float), np.nan).reshape(12, -1)[:, 0], 3)}')
