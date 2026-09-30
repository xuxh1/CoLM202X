#!/usr/bin/env python
"""Inundated wetland area by latitude band from two observational products and
the model: WAD2M v2.0 and GIEMS-MethaneCentric v1.1 (0.25 degree, 2010-2019),
and a V2 global run (wetland tile + floodplain, f_methane_finundated x active
area). Per band: mean of per-cell annual maxima and mean of monthly area (Mkm2),
the model also split into wetland tile and floodplain (soil tile). V2_MONTHS
(e.g. "6 7 8 9") restricts both to those months, since the satellites do not
see water under snow and ice. Usage: area_bands.py <ver> <name> <year>"""
import os
import sys

import netCDF4 as nc
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
BANDS = (('90S-30S', -90, -30), ('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90), ('all', -90, 90))
MONTHS = [int(m) - 1 for m in os.environ.get('V2_MONTHS', ' '.join(map(str, range(1, 13)))).split()]
PRODUCTS = {'WAD2M': ('/share/home/dq076/data/methane/WAD2M/WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc', 'Fw', 2000),
            'GIEMS-MC': ('/share/home/dq076/mode/Methane/data/GIEMS/GIEMS_MethaneCentric_1992-2020_025deg.nc', 'Fw', 1992)}


def product(path, var, y_first):
    with nc.Dataset(path) as h:
        lat = np.asarray(h['lat'][:])
        dl = abs(float(lat[1] - lat[0]))
        area = (np.radians(dl) ** 2) * 6371.0 ** 2 * np.cos(np.radians(lat))
        mx = {b: [] for b, *_ in BANDS}
        mn = {b: [] for b, *_ in BANDS}
        for y in range(2010, 2020):
            i0 = (y - y_first) * 12
            x = np.clip(np.ma.filled(h[var][i0:i0 + 12].astype(np.float64), 0.0), 0, 1)[MONTHS]
            x[~np.isfinite(x)] = 0.0
            pm, mm = x.max(axis=0), x.mean(axis=0)
            for b, lo, hi in BANDS:
                m = (lat >= lo) & (lat < hi)
                mx[b].append((pm[m] * area[m, None]).sum() / 1e6)
                mn[b].append((mm[m] * area[m, None]).sum() / 1e6)
    return {b: (np.mean(mx[b]), np.mean(mn[b])) for b in mx}


def model(ver, name, y):
    c = B.Case(ver, name)
    A_land, A_act, _ = VB.load_areas(c, y)
    with xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc') as tr:
        tr = tr.isel(time=MONTHS)
        fa = tr['f_methane_finundated'].fillna(0) * A_act
        fs = tr['f_methane_soil_finundated'].fillna(0) * tr['f_methane_area_soil'].fillna(0) * A_land
        out, parts = {}, {}
        for b, lo, hi in BANDS:
            m = (fa.lat >= lo) & (fa.lat < hi)
            s, sl = fa.where(m), fs.where(m)
            out[b] = (float(s.max('time').sum()) / 1e12, float(s.mean('time').sum()) / 1e12)
            parts[b] = (float((s - sl).mean('time').sum()) / 1e12, float(sl.mean('time').sum()) / 1e12)
    return out, parts


def main():
    ver, name, y = sys.argv[1], sys.argv[2], int(sys.argv[3])
    res = {k: product(*v) for k, v in PRODUCTS.items()}
    res[f'model {name} {y}'], parts = model(ver, name, y)
    print(f'# Mkm2: per-cell max (mean over years) / monthly mean over months {[m + 1 for m in MONTHS]}; '
          'products 2010-2019')
    print(f'{"band":8s} ' + ' | '.join(f'{k:>22s}' for k in res))
    for b, *_ in BANDS:
        print(f'{b:8s} ' + ' | '.join(f'{res[k][b][0]:10.2f} / {res[k][b][1]:9.2f}' for k in res)
              + f'   model monthly mean: wetland tile {parts[b][0]:.2f}, floodplain {parts[b][1]:.2f}')


if __name__ == '__main__':
    main()
