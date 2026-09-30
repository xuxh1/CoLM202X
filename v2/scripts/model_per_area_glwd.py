#!/usr/bin/env python
"""Wetland-tile CH4 emission per unit area by GLWD v2 wetland type, replacing
the five climate zones of get_wetland_veg_proxy that model_per_area.py splits
by (C-13 retires them). Each 2-degree cell goes to the type group holding most
of its GLWD tile area (v2/data/glwd_class_2deg.nc); marsh cells are also split
by band. Per group: tile area (Mkm2), emission (Tg yr-1), per wetland area and
per mean inundated area (g CH4 m-2 yr-1), mean inundated fraction.
Years from V2_YEARS (default 2003 2004). Usage: model_per_area_glwd.py <ver> <name>"""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
V2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
YEARS = tuple(int(y) for y in os.environ.get('V2_YEARS', '2003 2004').split())
GROUPS = {'swamp forest (16,18)': (16, 18), 'marsh (17,19)': (17, 19), 'boreal peat forested (22)': (22,),
          'boreal peat open (23)': (23,), 'temperate peat (24,25)': (24, 25),
          'tropical peat forested (26)': (26,), 'tropical peat open (27)': (27,)}


def dominant_group(lat, lon):
    """Group label per model cell (nearest GLWD cell; '' where no tile class)."""
    g = xr.open_dataset(f'{V2}/data/glwd_class_2deg.nc')
    a = np.stack([sum(g[f'area_class_{c}'].values for c in cls) for cls in GROUPS.values()])
    lab = np.array(list(GROUPS))[a.argmax(axis=0)]
    lab[a.sum(axis=0) <= 0] = ''
    out = xr.DataArray(lab, coords={'lat': g.lat, 'lon': g.lon}, dims=('lat', 'lon'))
    return out.reindex(lat=lat, lon=lon, method='nearest', tolerance=0.5).fillna('')


def main():
    ver, name = sys.argv[1], sys.argv[2]
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, YEARS[0])
    acc = {}
    for y in YEARS:
        with xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc') as tr:
            w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
            g = lambda v: tr[v].fillna(0)
            aw, asl = g('f_methane_area_wetland'), g('f_methane_area_soil')
            fw = ((g('f_methane_finundated') * (aw + asl) - g('f_methane_soil_finundated') * asl)
                  / aw.where(aw > 0)).clip(0, 1).fillna(0)
            d = {'Aw': (aw * A_land).mean('time'), 'FAw': (fw * aw * A_land).mean('time'),
                 'Ew': (g('f_methane_surf_flux_wetland') * A_land * w).sum('time') * B.M_CH4}
            for k, v in d.items():
                acc[k] = acc.get(k, 0) + v.load() / len(YEARS)
    lab = dominant_group(acc['Aw'].lat, acc['Aw'].lon)
    L = [f'== {ver}/{name}, years {YEARS} (g CH4 m-2 yr-1; A Mkm2; E Tg/yr); cell -> GLWD group with most tile area']
    rows = [(k, lab == k) for k in GROUPS]
    for band, lo, hi in (('30-60N', 30, 60), ('30S-30N', -30, 30)):
        rows.append((f'  marsh {band}', (lab == 'marsh (17,19)') & (lab.lat >= lo) & (lab.lat < hi)))
    # tropical peat by region (measured fluxes differ by region: Q-22)
    lon180 = ((lab.lon + 180) % 360) - 180
    trop = (lab == 'tropical peat forested (26)') | (lab == 'tropical peat open (27)')
    for reg, lo, hi in (('SE Asia', 90, 160), ('S America', -95, -30), ('Africa', -20, 55)):
        rows.append((f'  tropical peat {reg}', trop & (lon180 >= lo) & (lon180 < hi)))
    rows.append(('all', xr.ones_like(lab, dtype=bool)))
    for k, m in rows:
        m = m & (acc['Aw'] > 0)
        A = float(acc['Aw'].where(m).sum()) / 1e12
        if A <= 0:
            continue
        FA = float(acc['FAw'].where(m).sum()) / 1e12
        E = float(acc['Ew'].where(m).sum())
        L.append(f'  {k:28s} A {A:5.2f}  E {E:6.1f}  per wetland area {E / A:6.1f}  mean F {FA / A:4.2f}  '
                 f'per mean inundated area {E / FA if FA else np.nan:6.1f}')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
