#!/usr/bin/env python
"""C-48 input: monthly flood series for river-fed towers, from GIEMS-MC
(GIEMS_MethaneCentric_1992-2020_025deg.nc, monthly inundated fraction Fw of
each 0.25-degree cell). A tower in a floodplain sees its floodplain either
under water or not, while the cell fraction also counts the dry land around
it, so the series is the largest Fw of the 3x3 cells around the tower divided
by its largest value over 1992-2020 (the share of the floodplain under
water). Depth is taken as 0.5 m times that share (assumed; no depth data at
these towers). Writes v2/data/site_flood_giems.nc (dims site, time; lat,
lon, year, month, flood_frac, flood_depth), January 1992 to December 2020.
Usage: mk_site_flood.py SITE [SITE ...]"""
import glob
import os
import sys

import netCDF4 as nc
import numpy as np

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
G = '/share/home/dq076/data/methane/GIEMS/GIEMS_MethaneCentric_1992-2020_025deg.nc'
S = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Sitedata'
DEPTH_MAX = 0.5


def main():
    sites = sys.argv[1:]
    g = nc.Dataset(G)
    la, lo = g['lat'][:], g['lon'][:]
    nt = g['Fw'].shape[0]
    frac = np.zeros((len(sites), nt))
    slat, slon = [], []
    for k, s in enumerate(sites):
        d = nc.Dataset(glob.glob(f'{S}/{s}_*site.nc')[0])
        y, x = float(d['latitude'][:].ravel()[0]), float(d['longitude'][:].ravel()[0])
        slat.append(y)
        slon.append(x)
        i, j = int(np.argmin(abs(la - y))), int(np.argmin(abs(lo - x)))
        fw = np.ma.filled(g['Fw'][:, i - 1:i + 2, j - 1:j + 2].astype(float), np.nan).reshape(nt, -1)
        m = np.nanmax(fw, 1)
        frac[k] = np.nan_to_num(m / max(np.nanmax(m), 1e-6))
        clim = [frac[k][mo::12].mean() for mo in range(12)]
        print(f'{s} ({y:.2f},{x:.2f}) record max Fw {np.nanmax(m):.2f}; mean monthly share', ' '.join(f'{c:.2f}' for c in clim))
    out = f'{V2}/data/site_flood_giems.nc'
    with nc.Dataset(out, 'w') as o:
        o.createDimension('site', len(sites))
        o.createDimension('time', nt)
        o.createVariable('lat', 'f8', ('site',))[:] = slat
        o.createVariable('lon', 'f8', ('site',))[:] = slon
        o.createVariable('year', 'i4', ('time',))[:] = 1992 + np.arange(nt) // 12
        o.createVariable('month', 'i4', ('time',))[:] = 1 + np.arange(nt) % 12
        v = o.createVariable('flood_frac', 'f8', ('site', 'time'))
        v[:] = frac
        v.long_name = 'share of the floodplain under water: 3x3 max GIEMS-MC Fw over its 1992-2020 max'
        v = o.createVariable('flood_depth', 'f8', ('site', 'time'))
        v[:] = DEPTH_MAX * frac
        v.units = 'm'
        v.long_name = f'{DEPTH_MAX} m times flood_frac (assumed)'
        o.sites = ' '.join(sites)
        o.source = 'mk_site_flood.py (C-48, paper V2)'


if __name__ == '__main__':
    main()
