#!/usr/bin/env python
"""Why the soil-tile (floodplain + upland) CH4 emission rose from the paper
global run (g2_paper_fp2) to V-3 (g2_v2). Years 2000-2001 of both.
Per run, globally and by latitude band: soil-tile area; its flooded area
(soil-tile inundated fraction x soil-tile area) as annual mean and as the sum
of per-cell annual maxima; emission E, production P, oxidation O; E and P per
unit flooded area; E/P; river flood area from the routing output.
Usage: floodplain_compare.py"""
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
BANDS = (('all', -90, 90), ('90S-30S', -90, -30), ('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90))


def load(ver, name):
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, YEARS[0])
    acc = {}
    for y in YEARS:
        tr = xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        g = lambda v: tr[v].fillna(0)
        As = g('f_methane_area_soil') * A_land
        fl = g('f_methane_soil_finundated') * As
        d = {'A': As.mean('time'), 'Fmean': fl.mean('time'), 'Fmax': fl.max('time'),
             'E': (g('f_methane_surf_flux_soil') * A_land * w).sum('time') * B.M_CH4,
             'P': (g('f_methane_prod_tot_soil') * A_land * w).sum('time') * B.M_CH4,
             'O': (g('f_methane_oxid_tot_soil') * A_land * w).sum('time') * B.M_CH4}
        for k, v in d.items():
            acc[k] = acc.get(k, 0) + v.load() / len(YEARS)
        tr.close()
        u = xr.open_dataset(f'{c.dir}/history/{name}_hist_unitcat_{y}.nc')
        fa = u['f_floodarea'].fillna(0)
        acc['riv_mean'] = acc.get('riv_mean', 0) + float(fa.mean('time').sum()) / 1e12 / len(YEARS)
        acc['riv_max'] = acc.get('riv_max', 0) + float(fa.max('time').sum()) / 1e12 / len(YEARS)
        u.close()
    return acc


def main():
    L = ['# soil tile, 2000-2001 mean; areas Mkm2; E, P, O Tg CH4/yr; per flooded area (annual mean) g CH4 m-2 yr-1']
    for lab, ver, name in RUNS:
        d = load(ver, name)
        L.append(f'== {lab}   routing flood area: annual mean {d["riv_mean"]:.2f}, sum of per-cell max {d["riv_max"]:.2f} '
                 '(units as written by the routing output, /1e12)')
        for b, lo, hi in BANDS:
            m = (d['A'].lat >= lo) & (d['A'].lat < hi)
            s = lambda k: float(d[k].where(m).sum())
            A, Fm, Fx = s('A') / 1e12, s('Fmean') / 1e12, s('Fmax') / 1e12
            E, P, O = s('E'), s('P'), s('O')
            L.append(f'  {b:8s} soil A {A:6.2f}  flooded mean {Fm:5.2f} max-sum {Fx:5.2f}  E {E:6.1f}  P {P:6.1f}  '
                     f'O {O:5.1f}  E/Fmean {E / Fm if Fm else np.nan:6.1f}  P/Fmean {P / Fm if Fm else np.nan:6.1f}  '
                     f'E/P {E / P if P else np.nan:4.2f}')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
