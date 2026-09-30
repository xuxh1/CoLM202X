#!/usr/bin/env python
"""Per-GCP-region reference values that the budget workbook does not give,
written once to v2/data/region_refs.csv for region_eval.py:
  lake_johnson  Johnson et al. (2022) lakes only (diffusion + ebullition,
                ice-out, fall turnover), Tg CH4/yr: daily climatology x
                area of each 0.25-degree cell (the rates are per grid-cell
                area: x cell area gives 41.6 Tg/yr as in the paper; x lake
                area would give 12.9)
  wad2m_area    WAD2M v2.0 inundated area, 2010-2019 mean of the per-year sum
                of per-pixel annual maxima (Mkm2), the statistic of the
                GCP 7.9 Mkm2 (v2/目标与判据.md Q-1)
A 0.25-degree cell takes the region of the 1-degree mask cell it lies in;
cells on ocean mask cells go to the nearest land mask cell.
Usage: region_refs.py"""
import os

import netCDF4 as nc
import numpy as np

G = '/share/home/dq076/data/methane/GCP/GCP-CH4-2024/GCP2023_regions_mask.nc'
J = '/share/home/dq076/data/methane/Johnson2022_lake'
W = '/share/home/dq076/data/methane/WAD2M/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'region_refs.csv')
NREG = 18


def mask_on(lat, lon):
    """Region index (0-17) of each (lat, lon) point of a regular grid."""
    with nc.Dataset(G) as m:
        reg = m['regions'][:].astype(int)
        names = [str(s).strip().replace(' ', '_') for s in m['regions_name'][:]]
        mlat, mlon = m['lat'][:], m['lon'][:]
    mlon = np.where(mlon > 180, mlon - 360, mlon)
    # fill ocean mask cells with the nearest land region (index-space BFS)
    r = reg.copy()
    from scipy import ndimage
    ocean = r == NREG
    _, idx = ndimage.distance_transform_edt(ocean, return_indices=True)
    r = r[idx[0], idx[1]]
    ii = np.array([np.argmin(np.abs(mlat - a)) for a in lat])
    jj = np.array([np.argmin(np.abs(((mlon - (b - 360 if b > 180 else b)) + 180) % 360 - 180)) for b in lon])
    return names[:NREG], r[np.ix_(ii, jj)]


def lakes():
    with nc.Dataset(f'{J}/Lake_Area.nc') as h:
        lat, lon = h['lat'][:], h['lon'][:]
    R = 6371e3
    area = ((np.radians(0.25) ** 2) * R ** 2 * np.cos(np.radians(np.asarray(lat))))[:, None] \
        * np.ones((1, lon.size))
    tot = np.zeros_like(area)
    for f, v in (('Lake_CH4_Diff_Ebul_Emiss.nc', 'DiffusionEbullition_TotalLakes'),
                 ('Lake_CH4_Ice_out_Emiss.nc', 'IceOut_TotalLakes'),
                 ('Lake_CH4_Fall_Turnover_Emiss.nc', 'FallTurnover_TotalLakes')):
        with nc.Dataset(f'{J}/{f}') as h:
            x = h[v]
            for t0 in range(0, x.shape[0], 30):
                a = np.ma.filled(x[t0:t0 + 30].astype(np.float64), 0.0)
                a = np.nan_to_num(a); a[a < 0] = 0
                tot += a.sum(0)
    return lat, lon, tot * area / 1e12          # g m-2 day-1 x day x m2 -> Tg


def wad2m():
    with nc.Dataset(W) as h:
        lat = h['lat'][:] if 'lat' in h.variables else h['latitude'][:]
        lon = h['lon'][:] if 'lon' in h.variables else h['longitude'][:]
        fw = h['Fw']
        R = 6371.0
        dlat = abs(float(lat[1] - lat[0]))
        a2 = ((np.radians(dlat) ** 2) * R ** 2 * np.cos(np.radians(np.asarray(lat))))[:, None]
        acc = np.zeros((lat.size, lon.size))
        ny = 0
        for y in range(2010 - 2000, 2020 - 2000):
            x = np.clip(np.ma.filled(fw[y * 12:(y + 1) * 12].astype(np.float64), 0.0), 0, 1)
            acc += x.max(0) * a2 / 1e6
            ny += 1
    return lat, lon, acc / ny


def main():
    rows = {}
    lat, lon, lk = lakes()
    names, r = mask_on(lat, lon)
    print(f'Johnson 2022 lakes, global {lk.sum():.1f} Tg/yr')
    for k in range(NREG):
        rows.setdefault(names[k], {})['lake_johnson'] = lk[r == k].sum()
    lat, lon, wa = wad2m()
    names, r = mask_on(lat, lon)
    print(f'WAD2M 2010-2019 annual-max area, global {wa.sum():.2f} Mkm2')
    for k in range(NREG):
        rows[names[k]]['wad2m_area'] = wa[r == k].sum()
    with open(OUT, 'w') as f:
        f.write('region,lake_johnson,wad2m_area\n')
        for n in names[:NREG]:
            f.write(f'{n},{rows[n]["lake_johnson"]:.3f},{rows[n]["wad2m_area"]:.4f}\n')
    print(open(OUT).read())


if __name__ == '__main__':
    main()
