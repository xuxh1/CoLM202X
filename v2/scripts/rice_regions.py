#!/usr/bin/env python
"""Rice by GCP land region for a V2 global case: rainfed (PFT 61) and
irrigated (PFT 62) rice area from the case's landdata/diag: cropfrac_elm is
the share of each crop type within the cropland (TypeIndex = CFT number =
PFT - 14: rice 47, irrigated rice 48; sums to 1), times patchfrac_elm of the
cropland land type (IGBP 12); the rice CH4 emission and the season
from the tracer history, and the emission per rice area. The inventories give
the regional emission (region_eval.py); the areas are for comparison with
harvested-area statistics (FAO: about 160 Mha, IRRI: 55 % irrigated).
Usage: rice_regions.py <version dir under cases/> <case name> [y0 y1]"""
import os
import sys

import numpy as np
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
    cf = xr.open_dataset(f'{case.dir}/landdata/diag/cropfrac_elm_2005.nc')
    ti = list(cf['TypeIndex'].values)
    frac = cf['cropfrac_elm'].where(cf['cropfrac_elm'] > -1e30).fillna(0)
    pf = xr.open_dataset(f'{case.dir}/landdata/diag/patchfrac_elm_2005.nc')
    crop = pf['patchfrac_elm'].sel(TypeIndex=12).where(lambda x: x > -1e30).fillna(0).values
    cellA = A_land.values * crop
    a61 = frac.isel(TypeIndex=ti.index(47)).values * cellA / 1e10          # Mha
    a62 = frac.isel(TypeIndex=ti.index(48)).values * cellA / 1e10
    e, mon, amax = [], [], []
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        fl = tr['f_methane_surf_flux_rice'].fillna(0) * A_land
        e.append(((fl * w).sum('time') * B.M_CH4).values)
        mon.append((fl * w * B.M_CH4).values)
        amax.append((tr['f_methane_area_rice'].fillna(0).max('time') * A_land / 1e10).values)
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    E = np.mean(e, 0)
    M = np.mean(mon, 0)
    AR = np.mean(amax, 0)
    names, W = RE.region_weights(lat, lon)
    print(f'# {ver}/{name}, {years[0]}-{years[-1]}; area Mha, emission Tg CH4/yr, intensity g CH4 m-2 yr-1 of rice area')
    print(f'{"region":24s} {"rainfed":>7s} {"irrig":>6s} {"irr%":>4s} {"CH4 area":>8s} {"E":>6s} {"E/area":>6s}  months (share of annual)')
    tot = np.zeros(5)
    for r, n in enumerate(names):
        x = [np.nansum(W[r] * v) for v in (a61, a62, AR, E)]
        m = np.array([np.nansum(W[r] * M[k]) for k in range(M.shape[0])])
        if x[0] + x[1] < 0.05:
            continue
        tot += np.array(x + [0])
        share = ''.join(f'{v:4.0%}'.replace('%', '') for v in m / max(m.sum(), 1e-9))
        inten = x[3] / (x[2] * 1e10) * 1e12 if x[2] > 0 else np.nan
        print(f'{n:24s} {x[0]:7.1f} {x[1]:6.1f} {x[1] / (x[0] + x[1]):4.0%} {x[2]:8.1f} {x[3]:6.1f} {inten:6.1f}  {share}')
    print(f'{"total":24s} {tot[0]:7.1f} {tot[1]:6.1f} {tot[1] / (tot[0] + tot[1]):4.0%} {tot[2]:8.1f} {tot[3]:6.1f}')


if __name__ == '__main__':
    main()
