#!/usr/bin/env python
"""C-45 input: which 2-degree cells double-crop their irrigated rice, from
MIRCA2000 (Portmann et al. 2010, doi:10.1029/2008GB003435; Zenodo 7422506,
crop 3). Physical irrigated rice area = the largest of the 12 monthly growing
areas (5 arc-min, stored pixel-interleaved), harvested area = the annual
harvested-area grid; the cropping intensity is harvested / physical over the
2-degree cell. Cells are nearly all single or double: with rice_double = 1
where the intensity is at least 1.5, one season everywhere plus a second on
those cells gives 101.5 Mha of irrigated rice harvest against 103.1 in MIRCA
(scripts/_tmp/分析_260928_rice_double_share.py in the Methane repository).
Writes v2/data/rice_double_2deg.nc (lat 89..-89, lon -179..179, the model
grid): rice_double (0 or 1), intensity_irrigated, physical_irrigated_ha,
harvested_irrigated_ha.
Usage: mk_rice_double.py"""
import gzip
import os
import zipfile

import netCDF4 as nc
import numpy as np

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
M = '/share/home/dq076/data/methane/MIRCA2000'


def main():
    z = zipfile.ZipFile(f'{M}/monthly_growing_area_grids.zip')
    a = np.frombuffer(gzip.decompress(z.read('monthly_growing_area_grids/crop_03_irrigated_12.flt.gz')),
                      dtype='<f4').astype(np.float64)
    a = np.where(a < 0, 0, a)
    g = np.moveaxis(a.reshape(2160, 4320, 12), 2, 0)                 # month, lat (N->S), lon (W->E)
    phys = g.reshape(12, 90, 24, 180, 24).sum((2, 4)).max(0)          # ha per 2-degree cell
    harv = np.load(f'{M}/rice_30mn.npz')['irc'].reshape(90, 4, 180, 4).sum((1, 3))
    inten = np.where(phys > 0, harv / np.maximum(phys, 1e-9), 0.0)
    dbl = (inten >= 1.5).astype(np.float64)
    lat = 89.0 - 2.0 * np.arange(90)
    lon = -179.0 + 2.0 * np.arange(180)
    out = f'{V2}/data/rice_double_2deg.nc'
    with nc.Dataset(out, 'w') as d:
        d.createDimension('lat', 90)
        d.createDimension('lon', 180)
        d.createVariable('lat', 'f8', ('lat',))[:] = lat
        d.createVariable('lon', 'f8', ('lon',))[:] = lon
        for name, val, units, ln in (
                ('rice_double', dbl, '1', 'irrigated rice sown twice a year (MIRCA2000 intensity >= 1.5)'),
                ('intensity_irrigated', inten, '1', 'irrigated rice harvested / physical area'),
                ('physical_irrigated_ha', phys, 'ha', 'largest monthly irrigated rice growing area'),
                ('harvested_irrigated_ha', harv, 'ha', 'irrigated rice harvested area')):
            v = d.createVariable(name, 'f8', ('lat', 'lon'))
            v.units = units
            v.long_name = ln
            v[:] = val
        d.source = 'mk_rice_double.py (C-45, paper V2)'
    print(f'double cells {int(dbl.sum())}; physical irrigated {phys.sum()/1e6:.1f} Mha, in double cells '
          f'{phys[dbl > 0].sum()/1e6:.1f}; emulated harvest {(phys.sum() + phys[dbl > 0].sum())/1e6:.1f} '
          f'vs {harv.sum()/1e6:.1f} Mha')


if __name__ == '__main__':
    main()
