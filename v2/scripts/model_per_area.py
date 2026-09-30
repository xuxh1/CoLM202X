#!/usr/bin/env python
"""Model CH4 emission per unit area, permanent wetlands vs floodplain wetlands,
on the same bases as the literature: per wetland area (annual, including
non-flooded time), per mean inundated area, and for floodplains per maximum
flooded area (the floodplain's own extent). Years from V2_YEARS (default 2000 2001).
Permanent wetland = wetland tile, split by the model's wetland biome class
(f_methane_wetland_type: 1 tropical peat, 2 tropical floodplain, 3 temperate
marsh, 4 boreal fen, 5 boreal bog). Floodplain = soil tile in cells whose soil
tile floods, split by latitude band. Usage: model_per_area.py [ver name]..."""
import os
import sys

import numpy as np
import xarray as xr

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B

YEARS = tuple(int(y) for y in os.environ.get('V2_YEARS', '2000 2001').split())   # complete years only
CLS = {1: 'tropical peat', 2: 'tropical floodplain*', 3: 'temperate marsh', 4: 'boreal fen', 5: 'boreal bog'}
BANDS = (('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90), ('90S-30S', -90, -30))


def run(ver, name):
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, YEARS[0])
    acc = {}
    for y in YEARS:
        tr = xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        g = lambda v: tr[v].fillna(0)
        aw, asl = g('f_methane_area_wetland'), g('f_methane_area_soil')
        fw = ((g('f_methane_finundated') * (aw + asl) - g('f_methane_soil_finundated') * asl)
              / aw.where(aw > 0)).clip(0, 1).fillna(0)
        fs = g('f_methane_soil_finundated')
        d = {'Aw': (aw * A_land).mean('time'), 'FAw': (fw * aw * A_land).mean('time'),
             'Ew': (g('f_methane_surf_flux_wetland') * A_land * w).sum('time') * B.M_CH4,
             'Es': (g('f_methane_surf_flux_soil') * A_land * w).sum('time') * B.M_CH4,
             'FAs_mean': (fs * asl * A_land).mean('time'), 'FAs_max': (fs * asl * A_land).max('time'),
             'wtype': tr['f_methane_wetland_type'].max('time')}
        for k, v in d.items():
            acc[k] = acc.get(k, 0) + v.load() / len(YEARS)
        tr.close()
    L = [f'== {ver}/{name}  (g CH4 m-2 yr-1; areas Mkm2; E Tg/yr)',
         'permanent wetland (wetland tile) by model biome class:']
    t = acc['wtype'].round()
    for k, lab in CLS.items():
        m = (t == k) & (acc['Aw'] > 0)
        A = float(acc['Aw'].where(m).sum()) / 1e12
        if A <= 0:
            continue
        FA = float(acc['FAw'].where(m).sum()) / 1e12
        E = float(acc['Ew'].where(m).sum())
        L.append(f'  {lab:22s} A {A:5.2f}  E {E:6.1f}  per wetland area {E / A:6.1f}  mean F {FA / A:4.2f}  '
                 f'per mean inundated area {E / FA if FA else np.nan:6.1f}')
    A = float(acc['Aw'].sum()) / 1e12
    FA = float(acc['FAw'].sum()) / 1e12
    E = float(acc['Ew'].sum())
    L.append(f'  {"all":22s} A {A:5.2f}  E {E:6.1f}  per wetland area {E / A:6.1f}  mean F {FA / A:4.2f}  '
             f'per mean inundated area {E / FA:6.1f}')
    L.append('floodplain wetland (soil tile, cells that flood) by band:')
    fl = acc['FAs_max'] > 0
    for b, lo, hi in BANDS + (('all', -90, 90),):
        m = fl & (acc['Es'].lat >= lo) & (acc['Es'].lat < hi)
        Amax = float(acc['FAs_max'].where(m).sum()) / 1e12
        Amean = float(acc['FAs_mean'].where(m).sum()) / 1e12
        E = float(acc['Es'].where(m).sum())
        if Amax <= 0:
            continue
        L.append(f'  {b:8s} max flooded A {Amax:5.2f}  mean flooded A {Amean:5.2f}  E {E:6.1f}  '
                 f'per max flooded area {E / Amax:6.1f}  per mean flooded area {E / Amean:6.1f}')
    return '\n'.join(L)


def main():
    runs = list(zip(sys.argv[1::2], sys.argv[2::2])) or [('paper_v2/v260924c', 'g2_v2'), ('v260910c', 'g2_paper_fp2')]
    print('\n'.join(run(v, n) for v, n in runs))


if __name__ == '__main__':
    main()
