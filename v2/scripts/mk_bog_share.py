#!/usr/bin/env python
"""Candidate 32 global input: open-bog share of the wetland classes on the 2-degree
grid of v2/data/glwd_class_2deg_c38.nc, from BAWLD (Olefeldt et al. 2021, 0.5-degree
cell fractions of permafrost bog PEB, wet tundra WTU, marsh MAR, bog BOG, fen FEN;
data/_tmp/分析_260909_BAWLD对比/bawld_cells.csv). The bog share sets the deeper
water-table floor DEF_METHANE%wetland_max_wtd_bog_share of the model.
Share = BOG / (PEB + WTU + MAR + BOG + FEN), averaged over the 0.5-degree cells of
each 2-degree cell weighted by their wetland area, as mk_rainfed_share.py does for
the rain-fed (PEB) share. Outside the BAWLD domain the share is 0 (floor as before).
Writes v2/data/glwd_class_2deg_c38_bog.nc: a copy of glwd_class_2deg_c38.nc (all
its variables, rainfed_share included) plus bog_share. Refuses to overwrite.
Usage: mk_bog_share.py"""
import os
import shutil
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
SRC = f'{V2}/data/glwd_class_2deg_c38.nc'
DST = f'{V2}/data/glwd_class_2deg_c38_bog.nc'
BAWLD = '/share/home/dq076/mode/Methane/data/_tmp/分析_260909_BAWLD对比/bawld_cells.csv'
Z_CLASS, Z_BOG = 0.13, 0.30          # boreal class floor of g2_bf and default z_bog [m]


def wmean(x, w, m):
    """Weighted mean of x over mask m."""
    return float((x * w)[m].sum() / w[m].sum())


def main():
    if os.path.exists(DST):
        sys.exit(f'{DST} exists; remove it first')
    shutil.copyfile(SRC, DST)
    b = pd.read_csv(BAWLD)
    w5 = b[['b_peb', 'b_wtu', 'b_mar', 'b_bog', 'b_fen']].sum(axis=1)
    b = b[w5 > 0].copy()
    b['bs'] = b.b_bog / w5[w5 > 0]
    b['w'] = w5[w5 > 0] * np.cos(np.radians(b.lat))           # wetland area weight
    with nc.Dataset(DST, 'a') as d:
        lat, lon = d['lat'][:], d['lon'][:]
        dl = abs(float(lat[1] - lat[0])) / 2
        num = np.zeros((lat.size, lon.size))
        den = np.zeros_like(num)
        ii = np.array([np.argmin(np.abs(lat - x)) for x in b.lat])
        jj = np.array([np.argmin(np.abs(((lon - x) + 180) % 360 - 180)) for x in b.lon])
        ok = np.abs(lat[ii] - b.lat.values) <= dl
        np.add.at(num, (ii[ok], jj[ok]), (b.bs * b.w).values[ok])
        np.add.at(den, (ii[ok], jj[ok]), b.w.values[ok])
        bs = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        v = d.createVariable('bog_share', 'f8', ('lat', 'lon'))
        v.units = '1'
        v.long_name = 'open-bog share of wetland: BAWLD bog BOG/(PEB+WTU+MAR+BOG+FEN); 0 outside BAWLD'
        v[:] = bs

        rf = np.asarray(d['rainfed_share'][:], dtype=float)
        rf = np.where(np.isfinite(rf), rf, 0.0)
        at = np.asarray(d['area_tile_classes'][:], dtype=float)
        at = np.where(np.isfinite(at), at, 0.0)
    # the model's weight of the bog share within the fed share, and the implied floor
    fed = 1.0 - np.clip(rf, 0, 1)
    fb = np.minimum(np.clip(bs, 0, 1), fed)
    wb = np.where(fb > 0, fb / np.maximum(fed, 1e-12), 0.0)
    zf = (1 - wb) * Z_CLASS + wb * Z_BOG
    lat2 = np.repeat(np.asarray(lat)[:, None], lon.size, 1)
    print(f'BAWLD 0.5-degree cells, wetland-area weighted: BOG share {(b.bs * b.w).sum() / b.w.sum():.3f}, '
          f'PEB share {(b.b_peb / (b.b_peb + b.b_wtu + b.b_mar + b.b_bog + b.b_fen) * b.w).sum() / b.w.sum():.3f}')
    print(f'2-degree cells with a BAWLD share: {int((den > 0).sum())}; with bog_share > 0: {int((bs > 0).sum())}')
    print('tile-area-weighted means (area_tile_classes)          bog_share  rainfed_share  w_b    floor[m]')
    for name, m in [('global', at > 0), ('north of 45N', (at > 0) & (lat2 > 45)),
                    ('north of 50N', (at > 0) & (lat2 > 50)), ('BAWLD cells', (at > 0) & (den > 0))]:
        print(f'  {name:14s} tile area {at[m].sum() / 1e6:6.3f} Mkm2   {wmean(bs, at, m):.3f}      '
              f'{wmean(rf, at, m):.3f}          {wmean(wb, at, m):.3f}  {wmean(zf, at, m):.3f}')
    m = (at > 0) & (den > 0)
    q = np.percentile(bs[m], [25, 50, 75])
    print(f'  bog_share over BAWLD cells with tile area, unweighted quartiles: '
          f'{q[0]:.2f} {q[1]:.2f} {q[2]:.2f}; tile area with bog_share >= 0.4: '
          f'{at[m & (bs >= 0.4)].sum() / 1e6:.3f} Mkm2')
    print(f'  floor as if every cell had z_class {Z_CLASS} m and z_bog {Z_BOG} m (the class floor '
          f'differs outside the boreal classes)')
    print(f'wrote {DST}')


if __name__ == '__main__':
    main()
