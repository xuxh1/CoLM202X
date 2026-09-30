#!/usr/bin/env python
"""Lake CH4 budget of a V2 global case in a few GCP regions (Q-29): per year the
lake-area mean of production, oxidation, ebullition, surface diffusion, the
sediment-to-water flux (g CH4 m-2 yr-1) and the end-of-year sediment column
and water-column stocks (g CH4 m-2), to see whether stocks drift; then the
monthly climatology of ebullition and diffusion and of the sediment stock.
Usage: lake_budget.py <version dir under cases/> <case name> <y0> <y1> [region ...]"""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402
import region_eval as RE                                             # noqa: E402

B = VB.B
MOL = 16.043                    # g CH4 per mol


def main():
    ver, name, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    regs = sys.argv[5:] or ['Canada', 'Russia', 'Brazil', 'Equatorial_Africa']
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    _, _, A_lake = VB.load_areas(case, years[0])
    AL = A_lake.fillna(0).values
    rows, mon = {r: [] for r in regs}, {r: [] for r in regs}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        if y == years[0]:
            names, W = RE.region_weights(tr['lat'].values, tr['lon'].values)
        sec = B.SEC_PER_MONTH[:tr.sizes['time']]
        v = {k: tr[k].fillna(0).values for k in (
            'f_methane_prod_tot_lake', 'f_methane_oxid_tot_lake', 'f_methane_surf_ebul_lake',
            'f_methane_surf_diff_lake', 'f_lake_sed_ch4_flux', 'f_totcol_methane_lake',
            'f_lake_water_ch4_stock')}
        tr.close()
        for r in regs:
            w = W[names.index(r)] * AL
            aw = w.sum()

            def m(a):                   # monthly lake-area mean
                return (a * w[None]).sum((1, 2)) / aw
            fl = [np.sum(m(v[k]) * sec) * MOL for k in list(v)[:5]]
            st = [m(v[k])[-1] * MOL for k in list(v)[5:]]
            rows[r].append((y, *fl, *st))
            mon[r].append(np.array([m(v[k]) * MOL * 86400 * 30 for k in list(v)[2:4]] + [m(v['f_totcol_methane_lake']) * MOL]))
    for r in regs:
        print(f'# {r}: lake-area mean, fluxes g CH4 m-2 yr-1, stocks g CH4 m-2 at year end')
        print(f'{"year":>4s} {"P":>6s} {"O":>6s} {"ebul":>6s} {"diff":>6s} {"sed->w":>6s} {"sedstk":>7s} {"wstk":>6s}')
        for x in rows[r]:
            print(f'{x[0]:4d} ' + ' '.join(f'{v:6.2f}' for v in x[1:6]) + f' {x[6]:7.2f} {x[7]:6.3f}')
        c = np.mean(mon[r], 0)
        print('month       ' + ''.join(f'{i:6d}' for i in range(1, 13)))
        print('ebul/mon    ' + ''.join(f'{v:6.2f}' for v in c[0]))
        print('diff/mon    ' + ''.join(f'{v:6.2f}' for v in c[1]))
        print('sed stock   ' + ''.join(f'{v:6.1f}' for v in c[2]))


if __name__ == '__main__':
    main()
