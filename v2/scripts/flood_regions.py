#!/usr/bin/env python
"""Wetland tile and floodplain by GCP land region for a V2 global case: CH4
(Tg/yr) and area (Mkm2, per-cell annual maximum) of the permanent-wetland tile
and of the floodplain (soil tile in cells whose soil tile floods in the year,
its inundated area f_methane_soil_finundated x soil-tile area), against WAD2M
(annual-max area, region_refs.csv) and the 22 GCP model runs (median wetland
CH4, region_eval.py).
Usage: flood_regions.py <version dir under cases/> <case name> [y0 y1]"""
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402
import region_eval as RE                                             # noqa: E402

B = VB.B


def main():
    ver, name = sys.argv[1], sys.argv[2]
    y0, y1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2010, 2019)
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, _, _ = VB.load_areas(case, years[0])
    acc = {k: [] for k in ('et', 'ef', 'at', 'af')}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')

        def tg(v):
            return ((tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4).values

        fsoil = tr['f_methane_soil_finundated'].fillna(0)
        dry = (fsoil.max('time') < 1e-4).values
        acc['et'].append(tg('f_methane_surf_flux_wetland'))
        acc['ef'].append(np.where(dry, 0, tg('f_methane_surf_flux_soil')))
        acc['at'].append((tr['f_methane_area_wetland'].fillna(0).max('time') * A_land / 1e12).values)
        acc['af'].append(((fsoil * tr['f_methane_area_soil'].fillna(0)).max('time') * A_land / 1e12).values)
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    m = {k: np.mean(v, 0) for k, v in acc.items()}
    names, W = RE.region_weights(lat, lon)
    bu = RE.gcp_table('2010-2019 BU Budget', 'Wetlands')
    rf = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                                  'region_refs.csv')).set_index('region')
    print(f'# {ver}/{name}, {years[0]}-{years[-1]}; CH4 Tg/yr, area Mkm2 (annual max)')
    print(f'{"region":24s} {"E tile":>6s} {"E flood":>7s} {"E sum":>6s} {"BU med":>6s} | '
          f'{"A tile":>6s} {"A flood":>7s} {"WAD2M":>6s} | {"tile/A":>6s} {"flood/A":>7s}')
    for r, n in enumerate(names):
        x = {k: np.nansum(W[r] * v) for k, v in m.items()}
        print(f'{n:24s} {x["et"]:6.1f} {x["ef"]:7.1f} {x["et"] + x["ef"]:6.1f} {np.median(bu[n].dropna()):6.1f} | '
              f'{x["at"]:6.3f} {x["af"]:7.3f} {rf.loc[n, "wad2m_area"]:6.3f} | '
              f'{x["et"] / max(x["at"], 1e-9):6.1f} {x["ef"] / max(x["af"], 1e-9):7.1f}')


if __name__ == '__main__':
    main()
