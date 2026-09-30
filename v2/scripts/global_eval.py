#!/usr/bin/env python
"""Global check of a V2 case against the GCP targets (v2/目标与判据.md).
Emission per source from scripts/figures/budget/ch4_budget_core (land-area
weights; wetland = wetland tile, soil = soil tile incl. floodplain and upland
uptake), 2010-2019 mean, and the wetland-tile emission by latitude band.
Wetland = wetland tile + floodplain (the soil tile in cells whose soil tile floods;
the rest of the soil tile is reported as upland), v2/目标与判据.md 第 1 节.
Inundated area: f_methane_finundated x active area (wetland + soil tiles) and
the wetland tile alone (f_methane_finundated_wetland if written, else skipped);
per year the sum of per-cell annual maxima, as WAD2M / GCP report it.
Usage: global_eval.py <version dir under cases/> <case name> [y0 y1]"""
import glob
import os
import sys

import numpy as np
import xarray as xr

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B

BANDS = (('90S-30S', -90, -30), ('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90))


def main():
    ver, name = sys.argv[1], sys.argv[2]
    y0, y1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2010, 2019)
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, A_active, A_lake = VB.load_areas(case, years[0])
    rows = {y: B.year_budget(y, case, (A_land, A_active, A_lake)) for y in years}
    L = [f'# {ver}/{name}, {years[0]}-{years[-1]} ({len(years)} years); Tg CH4/yr']
    for k in ('E_wetland', 'E_soil', 'E_rice', 'E_lake', 'E_global'):
        v = np.array([rows[y][k] for y in years])
        L.append(f'{k:10s} mean {v.mean():7.1f}  min {v.min():7.1f}  max {v.max():7.1f}')
    band = {b: [] for b, *_ in BANDS}
    band_fp = {b: [] for b, *_ in BANDS}
    upland, amax_act, amax_wet = [], [], []
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        e = (tr['f_methane_surf_flux_wetland'] * A_land * w).sum('time') * B.M_CH4
        es = (tr['f_methane_surf_flux_soil'].fillna(0) * A_land * w).sum('time') * B.M_CH4
        # upland part of the soil tile: cells whose soil tile never floods in the year
        dry = tr['f_methane_soil_finundated'].fillna(0).max('time') < 1e-4
        upland.append(float(es.where(dry).sum()))
        for b, lo, hi in BANDS:
            m = (e.lat >= lo) & (e.lat < hi)
            band[b].append(float(e.where(m).sum()))
            band_fp[b].append(float(es.where(m & ~dry).sum()))
        fin = tr['f_methane_finundated'].fillna(0)
        amax_act.append(float((fin * A_active).max('time').sum()) / 1e12)
        if 'f_methane_area_wetland' in tr:
            aw = tr['f_methane_area_wetland'].fillna(0) * A_land
            amax_wet.append(float(aw.max('time').sum()) / 1e12)
        tr.close()
    fp = np.array([rows[y]['E_soil'] for y in years]) - np.array(upland)
    tot = np.array([rows[y]['E_wetland'] for y in years]) + fp
    L.append(f'E_floodplain (soil tile in cells that flood) mean {fp.mean():7.1f}; upland soil tile {np.mean(upland):.2f}')
    L.append(f'E_wetland_total = wetland tile + floodplain: mean {tot.mean():7.1f}  min {tot.min():7.1f}  max {tot.max():7.1f}')
    L.append('wetland-tile emission by band: ' + ', '.join(f'{b} {np.mean(v):.1f}' for b, v in band.items()))
    L.append('wetland tile + floodplain by band: ' + ', '.join(f'{b} {np.mean(band[b]) + np.mean(band_fp[b]):.1f}' for b in band))
    L.append(f'inundated area, sum of per-cell annual max of finundated x active area: '
             f'{np.mean(amax_act):.2f} Mkm2 (years {min(amax_act):.2f}-{max(amax_act):.2f})')
    if amax_wet:
        L.append(f'wetland tile area (static): {np.mean(amax_wet):.2f} Mkm2')
    L.append('targets (compare E_wetland_total): GCP median 161.2 (IQR 136.6-179.3, range 119-203); area 8.0 +- 2.0 Mkm2 '
             '(per-cell annual max); bands 90S-30S 3 [1-6], 30S-30N 110 [64-146], 30N-60N 34 [17-64], 60N-90N 12 [4-30]')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
