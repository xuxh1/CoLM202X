#!/usr/bin/env python
"""Rice CH4 budget by GCP region for a V2 global case: production, oxidation and
emission of the rice columns (Tg CH4/yr) and E/P, and the growing-season mean
of the patch-mean inundated fraction over rice cells (cells with a rice area
fraction above 1%), to tell whether low rice CH4 comes from dry paddies or
from little production.
Usage: rice_budget.py <version dir under cases/> <case name> <y0> <y1>"""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402
import region_eval as RE                                             # noqa: E402

B = VB.B


def main():
    ver, name, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, _, _ = VB.load_areas(case, years[0])
    acc = {k: [] for k in ('P', 'O', 'E', 'fin')}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')

        def tg(v):
            return ((tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4).values
        acc['P'].append(tg('f_methane_prod_tot_rice'))
        acc['O'].append(tg('f_methane_oxid_tot_rice'))
        acc['E'].append(tg('f_methane_surf_flux_rice'))
        # growing-season inundation: months where the rice flux is above 10% of its peak
        fl = tr['f_methane_prod_tot_rice'].fillna(0)
        season = fl > 0.1 * fl.max('time')
        fin = tr['f_methane_finundated'].fillna(0)
        acc['fin'].append(fin.where(season).mean('time').fillna(0).values)
        ar = tr['f_methane_area_rice'].fillna(0).max('time').values
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    m = {k: np.mean(v, 0) for k, v in acc.items()}
    names, W = RE.region_weights(lat, lon)
    rice = ar > 0.01
    print(f'# {ver}/{name} {years[0]}-{years[-1]}: rice CH4 Tg/yr; fin = growing-season patch-mean inundated fraction in rice cells (rice-area weighted)')
    print(f'{"region":24s} {"P":>6s} {"O":>6s} {"E":>6s} {"E/P":>5s} {"fin":>5s}')
    for r, n in enumerate(names):
        x = {k: np.nansum(W[r] * m[k]) for k in ('P', 'O', 'E')}
        if x['P'] < 0.05:
            continue
        wr = W[r] * ar * rice
        f = np.nansum(wr * m['fin']) / max(np.nansum(wr), 1e-12)
        print(f'{n:24s} {x["P"]:6.1f} {x["O"]:6.1f} {x["E"]:6.1f} {x["E"] / x["P"]:5.2f} {f:5.2f}')


if __name__ == '__main__':
    main()
