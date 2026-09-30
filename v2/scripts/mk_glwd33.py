#!/usr/bin/env python
"""All 34 GLWD v2 classes (0 dryland, 1-33 lakes, rivers and wetlands; Lehner et al.
2025) on the 2-degree model grid of v2/data/glwd_class_2deg_c38.nc and around the 77
flux towers of scripts/sites/LIST_sites_77.csv. Input for the continuous
class-weighted wetland typing that replaces the latitude / organic-matter split.

Source: GLWD_v2_0_area_by_class_pct, one GeoTIFF per class, uint8 percent of the
15-arcsecond pixel area (0-100), 255 = no data (ocean; coverage 84N-56S).

Output 1, v2/data/glwd33_2deg.nc: area_class_00 ... area_class_33 [km2] per 2-degree
cell = sum over its pixels of pct/100 x pixel area, the grid and convention of
glwd_class_2deg.py (values > 100 count as 0, spherical pixel area with R = 6371.0072
km, rows from 84N); area_land [km2] = area of the pixels with data in the class-00
layer (land including inland water). Layers 16-19 and 22-27 are checked against
glwd_class_2deg_c38.nc.

Output 2, v2/data/glwd33_sites77.csv: per tower (coordinates, IGBP and
SITE_CLASSIFICATION from the global attributes of
data/FLUXNET-CH4/Observation/<ID>_*_Flux.nc) the area share of each class in the disc
of 1 km and of 5 km radius around the tower (class area / land area in the disc; each
pixel weighted by its area times the part of it inside the disc, from a 10 x 10
sub-pixel sample), the land fraction of each disc, the dominant class with dryland
(main_*) and among classes 1-33 (main_wet_*, 0 if the disc has no wetland; ties go to
the lower class, as GLWD main_class), and the dominant class of the tower's own pixel
(pix_main, -1 if the pixel has no data).

Runs one worker per layer on a compute node, never on a login node. Refuses to
overwrite either output.
Usage: mk_glwd33.py [nproc]   (default 12 workers)"""
import glob
import multiprocessing as mp
import os
import socket
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
REF = f'{V2}/data/glwd_class_2deg_c38.nc'
DST_NC = f'{V2}/data/glwd33_2deg.nc'
DST_CSV = f'{V2}/data/glwd33_sites77.csv'
GLWD = '/share/home/dq076/mode/Methane/data/GLWD_v2/area_by_class_pct/'
TIF = GLWD + 'GLWD_v2_0_area_by_class_pct/GLWD_v2_0_class_{:02d}_pct.tif'
LEGEND = GLWD + 'GLWD_Legend_v2_0.csv'
SITES = '/share/home/dq076/mode/Methane/scripts/sites/LIST_sites_77.csv'
OBS = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Observation/'
FLXLIST = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/LIST_FLUXNET-CH4.csv'

K = list(range(34))
TILE = list(range(16, 20)) + list(range(22, 28))    # classes of the existing file
R = 6371.0072                                        # earth radius [km], as glwd_class_2deg.py
D = 1.0 / 240                                        # pixel size [deg]
NROW, NCOL, LAT_N = 33600, 86400, 84.0
CELL = 480                                           # pixels per 2-degree cell side
BLOCK = 1920                                         # rows per read: multiple of CELL and of the 128-row tile
RADII = (1.0, 5.0)                                   # disc radii around the towers [km]
NSUB = 10                                            # sub-pixel samples per pixel side
DISCS = []                                           # tower windows, set in main before the workers fork

# expected GLWD classes, used only for the printed consistency check of the towers
FOREST = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28}
PEAT, PALU = set(range(22, 28)), set(range(16, 20))
NONFOR = {9, 11, 13, 15, 17, 19, 21, 23, 25, 27}
TYPE = {'Bog': PEAT | PALU, 'Fen': PEAT | PALU, 'Wet tundra': PEAT | PALU | NONFOR,
        'Marsh': NONFOR | {30, 31}, 'Swamp': FOREST, 'Salt marsh': set(range(28, 33)),
        'Lake': {1, 2, 3, 6, 8, 9}, 'Rice': {33}}
DRAINED = {0, 33} | PEAT | PALU
IGBP_FOREST = {'ENF', 'EBF', 'DBF', 'DNF', 'MF'}
IGBP_OPEN = {'GRA', 'OSH', 'CSH', 'BSV', 'SNO', 'SAV'}


def row_area(i0, n):
    """Area [km2] of one pixel in each of the rows i0 .. i0+n-1."""
    lat_n = LAT_N - (i0 + np.arange(n)) * D
    return R ** 2 * np.radians(D) * (np.sin(np.radians(lat_n)) - np.sin(np.radians(lat_n - D)))


def read_sites():
    """The 77 towers with coordinates, IGBP and SITE_CLASSIFICATION from the Observation
    files, cross-checked against LIST_FLUXNET-CH4.csv and LIST_sites_77.csv."""
    s = pd.read_csv(SITES)
    fl = pd.read_csv(FLXLIST).set_index('ID')
    rows = []
    for sid in s.ID:
        f = glob.glob(f'{OBS}{sid}_*_Flux.nc')
        if len(f) != 1:
            sys.exit(f'{sid}: {len(f)} Observation files')
        with nc.Dataset(f[0]) as d:
            rows.append(dict(ID=sid, lat=float(d.getncattr('latitude')), lon=float(d.getncattr('longitude')),
                             IGBP=str(d.getncattr('IGBP_veg_short')),
                             SITE_CLASSIFICATION=str(d.getncattr('SITE_CLASSIFICATION'))))
    o = pd.DataFrame(rows)
    dd = np.maximum(np.abs(o.lat.values - fl.loc[o.ID, 'LAT'].values),
                    np.abs(o.lon.values - fl.loc[o.ID, 'LON'].values))
    print(f'towers {len(o)}; coordinates differing from LIST_FLUXNET-CH4.csv by > 1e-4 deg: '
          f'{list(o.ID[dd > 1e-4])}')
    m = (o.IGBP.values != s.IGBP.values) | (o.SITE_CLASSIFICATION.values != s.SITE_CLASSIFICATION.values)
    print(f'IGBP or SITE_CLASSIFICATION differing from LIST_sites_77.csv: {list(o.ID[m])}')
    return o


def site_discs(lat, lon):
    """Pixel window around each tower and, per radius, the weight [km2] of each pixel:
    pixel area x part of the pixel inside the disc (haversine distance of sub-pixel
    centres to the tower)."""
    out = []
    rmax, dy = max(RADII), R * np.radians(D)
    s = (np.arange(NSUB) + 0.5) / NSUB
    for la, lo in zip(lat, lon):
        ti, tj = int(np.floor((LAT_N - la) / D)), int(np.floor((lo + 180.0) / D))
        hr = int(np.ceil(rmax / dy)) + 2
        hc = int(np.ceil(rmax / (dy * np.cos(np.radians(min(abs(la) + hr * D, 89.0)))))) + 2
        r0, c0, h, w = ti - hr, tj - hc, 2 * hr + 1, 2 * hc + 1
        if r0 < 0 or c0 < 0 or r0 + h > NROW or c0 + w > NCOL:
            sys.exit(f'tower at {la}, {lo}: window leaves the raster')
        phi = np.radians(LAT_N - (r0 + np.arange(h)[:, None] + s[None, :]) * D)[:, :, None, None]
        lam = np.radians(-180.0 + (c0 + np.arange(w)[:, None] + s[None, :]) * D)[None, None, :, :]
        p0, l0 = np.radians(la), np.radians(lo)
        a = np.sin((phi - p0) / 2) ** 2 + np.cos(p0) * np.cos(phi) * np.sin((lam - l0) / 2) ** 2
        dist = 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))          # (h, NSUB, w, NSUB) [km]
        area = row_area(r0, h)[:, None]
        wts = [(dist <= r).mean(axis=(1, 3)) * area for r in RADII]
        out.append(dict(r0=r0, c0=c0, h=h, w=w, ti=ti - r0, tj=tj - c0, wts=wts))
    return out


def layer(k):
    """One pass over the layer of class k: per 2-degree cell the class area, the area
    and count of no-data pixels; per tower disc the class area and the data area; the
    value of the tower pixel."""
    nlat, nlon = NROW // CELL, NCOL // CELL
    area, valid = np.zeros((nlat, nlon)), np.zeros((nlat, nlon))
    nod, odd = np.zeros((nlat, nlon), dtype=np.int64), 0
    ns = len(DISCS)
    sa, sv = np.zeros((len(RADII), ns)), np.zeros((len(RADII), ns))
    pix = np.zeros(ns, dtype=np.int64)
    with rasterio.Env(GDAL_CACHEMAX=512), rasterio.open(TIF.format(k)) as src:
        if (src.width, src.height, src.dtypes[0], src.nodata) != (NCOL, NROW, 'uint8', 255.0):
            raise RuntimeError(f'class {k}: unexpected raster layout')
        t = src.transform
        if abs(t.a - D) > 1e-12 or abs(t.c + 180) > 1e-9 or abs(t.f - LAT_N) > 1e-9:
            raise RuntimeError(f'class {k}: unexpected transform {t}')
        for i0 in range(0, NROW, BLOCK):
            n = min(BLOCK, NROW - i0)
            p = src.read(1, window=Window(0, i0, NCOL, n))
            miss = p == 255
            bad = p > 100
            odd += int(bad.sum()) - int(miss.sum())                  # values 101-254, expected none
            p[bad] = 0
            ra = row_area(i0, n)[:, None]
            s = p.reshape(n, nlon, CELL).sum(axis=2, dtype=np.int64)            # pct sum per row and cell
            m = CELL - miss.reshape(n, nlon, CELL).sum(axis=2, dtype=np.int64)  # data pixels per row and cell
            c = slice(i0 // CELL, (i0 + n) // CELL)
            area[c] = (s * ra / 100.0).reshape(n // CELL, CELL, nlon).sum(axis=1)
            valid[c] = (m * ra).reshape(n // CELL, CELL, nlon).sum(axis=1)
            nod[c] = (CELL - m).reshape(n // CELL, CELL, nlon).sum(axis=1)
        for j, q in enumerate(DISCS):
            p = src.read(1, window=Window(q['c0'], q['r0'], q['w'], q['h'])).astype(np.float64)
            pix[j] = int(p[q['ti'], q['tj']])
            ok = p <= 100
            pp = np.where(ok, p, 0.0) / 100.0
            for r, wt in enumerate(q['wts']):
                sa[r, j], sv[r, j] = (wt * pp).sum(), (wt * ok).sum()
    print(f'  class {k:02d} read', flush=True)
    return k, area, valid, nod, odd, sa, sv, pix


def flags(row):
    """Disagreements of the 1-km dominant classes with SITE_CLASSIFICATION and IGBP."""
    out = []
    sc, ig, mn, mw = row.SITE_CLASSIFICATION, row.IGBP, row.main_1km, row.main_wet_1km
    if sc == 'Upland' and mn != 0:
        out.append(f'SITE Upland but class {mn} dominates')
    elif sc == 'Drained' and mn not in DRAINED:
        out.append(f'SITE Drained but class {mn} dominates')
    elif sc in TYPE:
        if mn == 0:
            out.append(f'SITE {sc} but dryland dominates')
        if mw not in TYPE[sc]:
            out.append(f'SITE {sc} but dominant wetland class {mw}')
    elif sc not in ('Upland', 'Drained'):
        out.append(f'SITE {sc} has no rule')
    if ig == 'WET' and mn == 0:
        out.append('IGBP WET but dryland dominates')
    elif ig == 'WAT' and not 1 <= mn <= 9:
        out.append(f'IGBP WAT but class {mn} dominates')
    elif ig == 'CRO' and mn not in (0, 33):
        out.append(f'IGBP CRO but class {mn} dominates')
    elif ig in IGBP_FOREST and mn > 0 and mn not in FOREST:
        out.append(f'IGBP {ig} but non-forested class {mn} dominates')
    elif ig in IGBP_OPEN and mn in FOREST:
        out.append(f'IGBP {ig} but forested class {mn} dominates')
    return '; '.join(out)


def main():
    global DISCS
    if socket.gethostname().startswith('mgt'):
        sys.exit('run on a compute node, not on a login node')
    for f in (DST_NC, DST_CSV):
        if os.path.exists(f):
            sys.exit(f'{f} exists; remove it first')
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    leg = pd.read_csv(LEGEND).set_index('GLWD_ID')['Class_name']
    st = read_sites()
    DISCS = site_discs(st.lat.values, st.lon.values)
    print(f'reading 34 layers with {nproc} workers', flush=True)
    with mp.get_context('fork').Pool(nproc) as pool:
        res = {r[0]: r[1:] for r in pool.imap_unordered(layer, K)}
    area = {k: res[k][0] for k in K}
    land, nod0 = res[0][1], res[0][2]

    # --- consistency of the layers ---------------------------------------------
    ndiff = {k: int((res[k][2] != nod0).sum()) for k in K}
    nodd = {k: res[k][3] for k in K}
    print(f'cells whose no-data pixel count differs from layer 00: '
          f'{ {k: v for k, v in ndiff.items() if v} or "none"}')
    print(f'pixels with values 101-254: {({k: v for k, v in nodd.items() if v}) or "none"}')
    tot = sum(area[k] for k in K)
    m = land > 1000
    r = tot[m] / land[m]
    print(f'land (pixels with data in layer 00) {land.sum() / 1e6:.3f} Mkm2; sum of the 34 class areas '
          f'{tot.sum() / 1e6:.3f} Mkm2 (ratio {tot.sum() / land.sum():.5f}); per cell with land > 1000 km2: '
          f'ratio min {r.min():.4f} max {r.max():.4f}')

    # --- grid and comparison with glwd_class_2deg_c38.nc ---------------------------
    with nc.Dataset(REF) as d:
        lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
        if not (np.allclose(lat, LAT_N - 1 - 2 * np.arange(NROW // CELL))
                and np.allclose(lon, -179 + 2 * np.arange(NCOL // CELL))):
            sys.exit('grid of the reference file differs from the 2-degree grid assumed here')
        print('comparison with glwd_class_2deg_c38.nc (km2)')
        print('  class   old Mkm2    new Mkm2    max|new-old| per cell   max|rel| (cells > 1 km2)')
        old_t = np.zeros_like(land)
        for k in TILE:
            old = np.asarray(d[f'area_class_{k:02d}'][:], float)
            old_t += old
            dif = np.abs(area[k] - old)
            big = old > 1
            rel = (dif[big] / old[big]).max() if big.any() else 0.0
            print(f'  {k:5d}   {old.sum() / 1e6:.6f}    {area[k].sum() / 1e6:.6f}    {dif.max():.3e}'
                  f'               {rel:.3e}')
        new_t = sum(area[k] for k in TILE)
        ref_t = np.asarray(d['area_tile_classes'][:], float)
        print(f'  16-19,22-27 total: old {ref_t.sum() / 1e6:.6f} Mkm2 (sum of old classes {old_t.sum() / 1e6:.6f}), '
              f'new {new_t.sum() / 1e6:.6f} Mkm2, max|new-old| per cell {np.abs(new_t - ref_t).max():.3e} km2')

    # --- global class table ------------------------------------------------------
    print('global area per class (Mkm2, % of land)')
    for k in K:
        print(f'  {k:2d}  {area[k].sum() / 1e6:7.3f}  {100 * area[k].sum() / land.sum():6.2f}  {leg[k]}')
    wet = sum(area[k] for k in K if k > 0)
    print(f'  classes 1-33: {wet.sum() / 1e6:.3f} Mkm2, {100 * wet.sum() / land.sum():.2f} % of land '
          f'(TechDoc: 18.2 Mkm2, 13.4 % of land excluding Antarctica)')

    # --- netCDF ----------------------------------------------------------------------
    with nc.Dataset(DST_NC, 'w') as o:
        o.createDimension('lat', lat.size)
        o.createDimension('lon', lon.size)
        v = o.createVariable('lat', 'f8', ('lat',))
        v[:], v.units, v.long_name = lat, 'degrees_north', 'latitude of the 2-degree cell centre'
        v = o.createVariable('lon', 'f8', ('lon',))
        v[:], v.units, v.long_name = lon, 'degrees_east', 'longitude of the 2-degree cell centre'
        for k in K:
            v = o.createVariable(f'area_class_{k:02d}', 'f8', ('lat', 'lon'))
            v[:], v.units, v.long_name = area[k], 'km2', f'area of GLWD v2 class {k}: {leg[k]}'
        v = o.createVariable('area_land', 'f8', ('lat', 'lon'))
        v[:], v.units = land, 'km2'
        v.long_name = 'area of the GLWD v2 pixels with data (land including inland water; ocean excluded)'
        o.source = ('GLWD v2.0 (Lehner et al. 2025), GLWD_v2_0_area_by_class_pct: per-class percent of the '
                    '15-arcsecond pixel area, 84N-56S; data DOI 10.6084/m9.figshare.28519994 (GLWD TechDoc v2.0)')
        o.method = ('area_class_kk = sum over the 480 x 480 pixels of each 2-degree cell of pct/100 x pixel '
                    'area; values > 100 (no data 255) count as 0; spherical pixel area with R = 6371.0072 km; '
                    'grid and convention of glwd_class_2deg.py (glwd_class_2deg_c38.nc). area_land = area of '
                    'the pixels with data (value <= 100) in the class-00 layer.')
        o.script = 'v2/scripts/mk_glwd33.py'

    # --- towers ------------------------------------------------------------------------
    ns = len(st)
    pv = np.array([res[k][6] for k in K], dtype=np.int64)            # (34, ns) tower-pixel percent
    pok = pv[0] != 255
    pvc = np.where(pv <= 100, pv, -1)
    cols = {c: st[c].values for c in ['ID', 'lat', 'lon', 'IGBP', 'SITE_CLASSIFICATION']}
    cols['pix_main'] = np.where(pok, pvc.argmax(0), -1)
    cols['pix_main_wet'] = np.where(pvc[1:].max(0) > 0, pvc[1:].argmax(0) + 1, 0)
    cols['pix_pct_main'] = np.where(pok, pvc.max(0), -1)
    print(f'tower pixels without data: {list(st.ID[~pok])}; percent sum of the tower pixels: '
          f'min {pvc.clip(0).sum(0)[pok].min()} max {pvc.clip(0).sum(0)[pok].max()}')
    for r, rad in enumerate(RADII):
        lab = f'{int(rad)}km'
        cls = np.array([res[k][4][r] for k in K])                      # (34, ns) class area in disc
        lnd = res[0][5][r]                                             # data area in disc
        dsk = np.array([q['wts'][r].sum() for q in DISCS])             # disc area
        sh = np.where(lnd > 0, cls / np.where(lnd > 0, lnd, 1.0), np.nan)
        sh0 = np.nan_to_num(sh)
        cols[f'land_frac_{lab}'] = np.round(lnd / dsk, 4)
        cols[f'main_{lab}'] = np.where(lnd > 0, sh0.argmax(0), -1)
        cols[f'main_wet_{lab}'] = np.where(sh0[1:].max(0) > 0, sh0[1:].argmax(0) + 1, 0)
        for k in K:
            cols[f'share_{lab}_{k:02d}'] = np.round(sh[k], 6)
        ss = sh0.sum(0)[lnd > 0]
        print(f'{lab} discs: area {dsk.min():.3f}-{dsk.max():.3f} km2 (exact {np.pi * rad ** 2:.3f}); '
              f'sum of the 34 shares min {ss.min():.4f} max {ss.max():.4f}')
    out = pd.DataFrame(cols)
    out.to_csv(DST_CSV, index=False)

    print('towers per dominant class within 1 km (with dryland | among classes 1-33)')
    c1, c2 = out.main_1km.value_counts(), out.main_wet_1km.value_counts()
    for k in sorted(set(c1.index) | set(c2.index)):
        print(f'  {k:2d}  {c1.get(k, 0):3d} | {c2.get(k, 0):3d}  {leg.get(k, "no data")}')
    out['wet_1km'] = 1 - out.share_1km_00
    out['flag'] = out.apply(flags, axis=1)
    print('towers whose 1-km dominant class disagrees with SITE_CLASSIFICATION or IGBP')
    print('  ID       SITE         IGBP  main  main_wet  wet share  pix  main_5km  reason')
    for _, q in out[out.flag != ''].iterrows():
        print(f'  {q.ID:8s} {q.SITE_CLASSIFICATION:12s} {q.IGBP:4s}  {q.main_1km:4d}  {q.main_wet_1km:8d}'
              f'  {q.wet_1km:9.2f}  {q.pix_main:3d}  {q.main_5km:8d}  {q.flag}')
    print(f'  {int((out.flag != "").sum())} of {ns} towers flagged')
    print(f'wrote {DST_NC}\nwrote {DST_CSV}')


if __name__ == '__main__':
    main()
