#!/usr/bin/env python
"""Monthly CH4 emission of a V2 global run by latitude band (Tg CH4 per month;
wetland tile + floodplain, as global_eval.py), with the share of the year in
the cold season (October-May) and in March-May.
Usage: band_monthly.py <ver> <name> <year>"""
import os
import sys

import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
BANDS = (('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90))


def main():
    ver, name, y = sys.argv[1], sys.argv[2], int(sys.argv[3])
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, y)
    with xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc') as tr:
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        g = lambda v: tr[v].fillna(0)
        E = (g('f_methane_surf_flux_wetland') + g('f_methane_surf_flux_soil')) * A_land * w * B.M_CH4   # Tg, as model_per_area.py
        print(f'# {ver}/{name} {y}: Tg CH4 per month, wetland tile + floodplain')
        print('band     ' + ' '.join(f'{m:6d}' for m in range(1, 13)) + '   year  Oct-May  Mar-May')
        for b, lo, hi in BANDS:
            s = E.where((E.lat >= lo) & (E.lat < hi)).sum(('lat', 'lon')).values
            tot = s.sum()
            cold = s[[9, 10, 11, 0, 1, 2, 3, 4]].sum() / tot
            spring = s[[2, 3, 4]].sum() / tot
            print(f'{b:8s} ' + ' '.join(f'{v:6.2f}' for v in s) + f'  {tot:5.1f}  {cold:7.2f}  {spring:7.2f}')


if __name__ == '__main__':
    main()
