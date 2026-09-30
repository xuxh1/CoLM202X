#!/usr/bin/env python
"""C-39: rice irrigated / rainfed split of the CFT surface data from MIRCA2000.
The 0.25-degree PCT_CFT of global_CFT_surface_data.nc (share of the crop land
unit per CFT) keeps its total rice share p47 + p48 (CFT 47 rainfed rice =
PFT 61, CFT 48 irrigated rice = PFT 62) in every cell; the irrigated part
becomes that total times the MIRCA2000 irrigated share of rice area of the
0.5-degree cell holding it (Portmann et al. 2010, doi:10.1029/2008GB003435;
Zenodo 7422506, crop 3). Cells where MIRCA has no rice keep the original
split. Writes the corrected file and a rawdata directory of symlinks to the
original with only this file replaced.
share = harvested (C-39): irrigated / total harvested area, which counts a
  double-cropped irrigated field twice (irrigated 62% globally);
share = physical (C-39b, default): irrigated / total physical area, the largest
  of the 12 monthly growing areas (irrigated 51% globally), for use with the
  second rice season of C-45 that adds the second harvest itself.
Usage: mk_cft_mirca_rice.py <original rawdata dir> <new rawdata dir> [physical|harvested]"""
import gzip
import os
import sys
import zipfile

import netCDF4 as nc
import numpy as np


def main():
    src, dst = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/')
    kind = sys.argv[3] if len(sys.argv) > 3 else 'physical'
    os.makedirs(dst, exist_ok=False)
    for e in sorted(os.listdir(src)):
        if e != 'global_CFT_surface_data.nc':
            os.symlink(f'{src}/{e}', f'{dst}/{e}')
    m = np.load('/share/home/dq076/data/methane/MIRCA2000/rice_30mn.npz')
    irc, rfc, mlat, mlon = m['irc'], m['rfc'], m['lat'], m['lon']
    if kind == 'physical':
        z = zipfile.ZipFile('/share/home/dq076/data/methane/MIRCA2000/monthly_growing_area_grids.zip')

        def phys(name):
            a = np.frombuffer(gzip.decompress(z.read(f'monthly_growing_area_grids/{name}')),
                              dtype='<f4').astype(np.float64)
            a = np.where(a < 0, 0, a)                      # 5 arc-min, 12 months pixel-interleaved
            return np.moveaxis(a.reshape(2160, 4320, 12), 2, 0).reshape(12, 360, 6, 720, 6).sum((2, 4)).max(0)
        irc, rfc = phys('crop_03_irrigated_12.flt.gz'), phys('crop_03_rainfed_12.flt.gz')
    tot = irc + rfc
    share = np.where(tot > 0, irc / np.maximum(tot, 1e-12), np.nan)
    fin, fout = f'{src}/global_CFT_surface_data.nc', f'{dst}/global_CFT_surface_data.nc'
    with nc.Dataset(fin) as a, nc.Dataset(fout, 'w', format=a.data_model) as b:
        b.setncatts(a.__dict__)
        b.setncattr('history_v2', f'C-39: CFT 47/48 rice split from MIRCA2000 irrigated share of {kind} area (mk_cft_mirca_rice.py)')
        for n, d in a.dimensions.items():
            b.createDimension(n, None if d.isunlimited() else len(d))
        for n, v in a.variables.items():
            w = b.createVariable(n, v.datatype, v.dimensions,
                                 fill_value=getattr(v, '_FillValue', None))
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != '_FillValue'})
            if n != 'PCT_CFT':
                w[:] = v[:]
        lat, lon = a['lat'][:], a['lon'][:]
        ii = np.array([np.argmin(np.abs(mlat - x)) for x in lat])
        jj = np.array([np.argmin(np.abs(((mlon - x) + 180) % 360 - 180)) for x in lon])
        s = share[np.ix_(ii, jj)]
        p = a['PCT_CFT'][:]
        p47 = np.ma.filled(p[:, :, 46].astype(np.float64), 0.0)
        p48 = np.ma.filled(p[:, :, 47].astype(np.float64), 0.0)
        r = p47 + p48
        use = (r > 0) & np.isfinite(s)
        n48 = np.where(use, r * np.nan_to_num(s), p48)
        n47 = np.where(use, r - n48, p47)
        q = np.ma.filled(p.astype(np.float64), 0.0)
        q[:, :, 46], q[:, :, 47] = n47, n48
        b['PCT_CFT'][:] = q
        w = np.cos(np.radians(np.asarray(lat)))[:, None]
        print(f'irrigated rice share of PCT_CFT (cos-lat weighted, no crop fraction): '
              f'{(p48 * w).sum() / (r * w).sum():.3f} -> {(n48 * w).sum() / (r * w).sum():.3f}; '
              f'cells changed {int(use.sum())}; max |sum change| {np.abs(q.sum(2) - np.ma.filled(p.astype(np.float64), 0.0).sum(2)).max():.2e}')


if __name__ == '__main__':
    main()
