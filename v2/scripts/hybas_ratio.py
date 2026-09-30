"""Local contributing-area ratio of the wetland around each flux tower from the
HydroBASINS level-12 sub-basin that holds the tower (D-34, preparing option 2
of Q-25; the GLWD-patch version of Q-18 failed where towers had no mapped
wetland nearby or sat in peatland complexes of thousands of km2).
Per tower: sub-basin area and upstream area (HydroBASINS SUB_AREA, UP_AREA),
and inside the sub-basin the GLWD v2 areas (class percent x pixel area) of the
model's wetland classes (16-19, 22-27), other wetland classes, open water
(1-7) and dryland (0); ratio R = dryland / model wetland, the upland that can
feed the wetland per unit wetland area. The required ratio of the earlier
offline water balance (v2/results/step2_ratio_check.log, tau30 column) is
printed beside it where it exists.
Usage: hybas_ratio.py"""
import csv
import glob
import re

import numpy as np
import pyogrio
import rasterio
import rasterio.features
import rasterio.windows
from shapely.geometry import Point

M = '/share/home/dq076/mode/Methane'
HB = '/share/home/dq076/data/HydroBASINS/v1c'
GL = f'{M}/data/GLWD_v2/area_by_class_pct/GLWD_v2_0_area_by_class_pct/GLWD_v2_0_class_{{:02d}}_pct.tif'
V2 = f'{M}/CoLM202X-paper-v2/v2'
WET = [16, 17, 18, 19, 22, 23, 24, 25, 26, 27]
WATER = [1, 2, 3, 4, 5, 6, 7]
OTHER = [8, 9, 10, 11, 12, 13, 14, 15, 20, 21, 28, 29, 30, 31, 32, 33]
REGIONS = ['eu', 'na', 'si', 'ar', 'as', 'af', 'sa', 'au']


def towers():
    ids = []
    for f in (f'{V2}/config/v2_k4lrxmnqf/LIST_sites.csv', f'{V2}/config/LIST_sites_val23.csv'):
        for r in csv.DictReader(open(f, newline='')):
            ids.append((r.get('ID') or r.get('site') or list(r.values())[0]).strip())
    out = []
    import netCDF4 as nc
    for s in dict.fromkeys(ids):
        p = glob.glob(f'{M}/data/FLUXNET-CH4/Sitedata/{s}_*_site.nc')
        if not p:
            continue
        with nc.Dataset(p[0]) as d:
            out.append((s, float(np.ravel(d['latitude'][:])[0]), float(np.ravel(d['longitude'][:])[0])))
    return out


def needed():
    d = {}
    for line in open(f'{V2}/results/step2_ratio_check.log'):
        m = re.match(r'^(\S{2}-\S{3})\s+\S+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+(\S+)\s+(\S+)\s*$', line)
        if m:
            d[m.group(1)] = float(m.group(3))
    return d


def basin(lat, lon):
    pt = Point(lon, lat)
    for r in REGIONS:
        path = f'zip://{HB}/hybas_{r}_lev12_v1c.zip!hybas_{r}_lev12_v1c.shp'
        gdf = pyogrio.read_dataframe(path, bbox=(lon - 0.01, lat - 0.01, lon + 0.01, lat + 0.01))
        hit = gdf[gdf.contains(pt)]
        if len(hit):
            return hit.iloc[0]
    return None


def areas(geom, lat):
    out = {}
    with rasterio.open(GL.format(0)) as src:
        win = rasterio.windows.from_bounds(*geom.bounds, transform=src.transform).round_offsets().round_lengths()
        win = win.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
        tr = src.window_transform(win)
        mask = rasterio.features.geometry_mask([geom], out_shape=(int(win.height), int(win.width)),
                                               transform=tr, invert=True, all_touched=False)
        rows = np.arange(int(win.height)) + 0.5
        lats = tr.f + rows * tr.e
        px = (abs(tr.a) * 111.32) * (abs(tr.e) * 111.32) * np.cos(np.radians(lats))[:, None]   # km2
    for group, cls in (('wet', WET), ('water', WATER), ('other', OTHER), ('dry', [0])):
        tot = 0.
        for c in cls:
            with rasterio.open(GL.format(c)) as src:
                a = src.read(1, window=win).astype(float)
                a[(a < 0) | (a > 100)] = 0.
            tot += float(np.sum(a / 100. * px * mask))
        out[group] = tot
    return out


def main():
    need = needed()
    print('site     lat     lon | HYBAS_ID      SUB_AREA  UP_AREA km2 | wet   other water  dry km2 | R=dry/wet | needed(tau30)')
    for s, lat, lon in towers():
        b = basin(lat, lon)
        if b is None:
            print(f'{s:7s} {lat:6.2f} {lon:7.2f} | no level-12 basin')
            continue
        a = areas(b.geometry, lat)
        R = a['dry'] / a['wet'] if a['wet'] > 0 else np.inf
        print(f"{s:7s} {lat:6.2f} {lon:7.2f} | {int(b['HYBAS_ID']):12d} {b['SUB_AREA']:8.1f} {b['UP_AREA']:9.1f} | "
              f"{a['wet']:6.1f} {a['other']:6.1f} {a['water']:5.1f} {a['dry']:6.1f} | {R:8.2f} | {need.get(s, np.nan):7.2f}")


if __name__ == '__main__':
    main()
