#!/usr/bin/env python
"""Why the wetland-tile emission per unit wetland area fell from the paper
global run (g2_paper_fp2: old surface data, saturated wetland tile, biome f)
to V-3 (g2_v2: GLWD surface data, dynamic wetland, f = 0.2, C-8, C-11).
Years 2000-2001 of both. Per run: wetland area, emission E, production P,
oxidation, E/A, P/A, E/P, and the wetland tile's inundated fraction
(active-area inundated fraction minus the soil tile's share), globally, by
latitude band, and split into cells that were wetland in the paper run
('old cells') and cells that are wetland only in the GLWD data ('new cells').
Usage: intensity_compare.py"""
import os
import sys

import numpy as np
import xarray as xr

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B

RUNS = (('paper', 'v260910c', 'g2_paper_fp2'), ('V-3', 'paper_v2/v260924c', 'g2_v2'))
YEARS = (2000, 2001)
BANDS = (('90S-30S', -90, -30), ('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90))


def load(ver, name):
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, YEARS[0])
    acc = {}
    for y in YEARS:
        tr = xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        g = lambda v: tr[v].fillna(0)
        aw, asl = g('f_methane_area_wetland'), g('f_methane_area_soil')
        # wetland-tile inundated fraction: active-area fraction minus the soil tile's part
        fw = ((g('f_methane_finundated') * (aw + asl) - g('f_methane_soil_finundated') * asl)
              / aw.where(aw > 0)).clip(0, 1)
        d = {'A': (aw * A_land).mean('time'),
             'E': (g('f_methane_surf_flux_wetland') * A_land * w).sum('time') * B.M_CH4,
             'P': (g('f_methane_prod_tot_wetland') * A_land * w).sum('time') * B.M_CH4,
             'O': (g('f_methane_oxid_tot_wetland') * A_land * w).sum('time') * B.M_CH4,
             'FA': (fw.fillna(0) * aw * A_land).mean('time')}
        for k, v in d.items():
            acc[k] = acc.get(k, 0) + v.load() / len(YEARS)
        tr.close()
    return acc


def line(lab, d, m):
    A = float(d['A'].where(m).sum()) / 1e12                          # Mkm2
    s = lambda k: float(d[k].where(m).sum())                         # Tg CH4 / yr (M_CH4 converts to Tg)
    E, P, O = s('E'), s('P'), s('O')                                 # Tg / Mkm2 = g m-2
    F = float(d['FA'].where(m).sum()) / 1e12 / A if A > 0 else np.nan
    return (f'{lab:22s} A {A:5.2f} Mkm2  E {E:6.1f}  P {P:6.1f}  O {O:5.1f}  E/A {E / A if A else np.nan:5.1f}  '
            f'P/A {P / A if A else np.nan:5.1f} g/m2/yr  E/P {E / P if P else np.nan:4.2f}  F_wet {F:4.2f}')


def main():
    D = {lab: load(ver, name) for lab, ver, name in RUNS}
    old = D['paper']['A'] > 0
    L = ['# wetland tile, 2000-2001 mean; E, P, O in Tg CH4/yr; E/A, P/A in g CH4 m-2 yr-1; '
         'F_wet = area-weighted mean inundated fraction of the wetland tile']
    for lab, _, _ in RUNS:
        d = D[lab]
        L.append(f'== {lab}')
        L.append(line('all', d, d['A'] > 0))
        if lab != 'paper':
            L.append(line('old cells', d, (d['A'] > 0) & old))
            L.append(line('new cells (GLWD only)', d, (d['A'] > 0) & ~old))
        for b, lo, hi in BANDS:
            m = (d['A'] > 0) & (d['A'].lat >= lo) & (d['A'].lat < hi)
            L.append(line(f'  {b}', d, m))
    print('\n'.join(L))


if __name__ == '__main__':
    main()
