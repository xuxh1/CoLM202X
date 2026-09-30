#!/usr/bin/env python
"""C-38 global input: rain-fed (ombrotrophic) share of the wetland classes on the
2-degree grid of v2/data/glwd_class_2deg.nc, from BAWLD (Olefeldt et al. 2021,
0.5-degree cell fractions of permafrost bog PEB, wet tundra WTU, marsh MAR, bog
BOG, fen FEN; data/_tmp/分析_260909_BAWLD对比/bawld_cells.csv). D-40: only the
permafrost bogs (palsas, peat plateaus, raised by ground ice above the
surrounding water table) count as cut off from lateral inflow; open bogs keep
the floor, which stands in for the peat water retention the model lacks (V-70).
Share = PEB / (PEB + WTU + MAR + BOG + FEN), averaged over the 0.5-degree cells
of each 2-degree cell weighted by their wetland area. Outside the BAWLD domain
the share is 0 (all wetland keeps its lateral inflow, as before). Writes
v2/data/glwd_class_2deg_c38.nc: a copy of glwd_class_2deg.nc plus rainfed_share.
Usage: mk_rainfed_share.py"""
import os
import shutil

import netCDF4 as nc
import numpy as np
import pandas as pd

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
SRC = f'{V2}/data/glwd_class_2deg.nc'
DST = f'{V2}/data/glwd_class_2deg_c38.nc'
BAWLD = '/share/home/dq076/mode/Methane/data/_tmp/分析_260909_BAWLD对比/bawld_cells.csv'


def main():
    shutil.copyfile(SRC, DST)
    b = pd.read_csv(BAWLD)
    w5 = b[['b_peb', 'b_wtu', 'b_mar', 'b_bog', 'b_fen']].sum(1)
    b = b[w5 > 0].copy()
    b['rf'] = b.b_peb / w5[w5 > 0]
    b['w'] = w5[w5 > 0] * np.cos(np.radians(b.lat))           # wetland area weight
    with nc.Dataset(DST, 'a') as d:
        lat, lon = d['lat'][:], d['lon'][:]
        dl = abs(float(lat[1] - lat[0])) / 2
        num = np.zeros((lat.size, lon.size))
        den = np.zeros_like(num)
        ii = np.array([np.argmin(np.abs(lat - x)) for x in b.lat])
        jj = np.array([np.argmin(np.abs(((lon - x) + 180) % 360 - 180)) for x in b.lon])
        ok = np.abs(lat[ii] - b.lat.values) <= dl
        np.add.at(num, (ii[ok], jj[ok]), (b.rf * b.w).values[ok])
        np.add.at(den, (ii[ok], jj[ok]), b.w.values[ok])
        rf = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        v = d.createVariable('rainfed_share', 'f8', ('lat', 'lon'))
        v.units = '1'
        v.long_name = 'share of wetland without lateral inflow: BAWLD permafrost bog PEB/(PEB+WTU+MAR+BOG+FEN); 0 outside BAWLD'
        v[:] = rf
        at = d['area_tile_classes'][:]
        at = np.where(np.isfinite(at), at, 0)
        print(f'cells with a share: {int((den > 0).sum())}; tile-area-weighted share north of 45N: '
              f'{(rf * at)[lat > 45].sum() / at[lat > 45].sum():.2f}')


if __name__ == '__main__':
    main()
