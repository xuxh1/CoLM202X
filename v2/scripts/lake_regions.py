#!/usr/bin/env python
"""Lakes by GCP land region for a V2 global case against Johnson et al. (2022):
model lake area (area_lake of the history, the land-water-body patches) and
lake CH4 and its production, oxidation, ebullition and diffusion per lake area
(set C variables, g CH4 m-2 yr-1), Johnson lake area (Lake_Area.nc, HydroLAKES-based) and lake CH4
(region_refs.csv), and the emission per lake area of both.
Usage: lake_regions.py <version dir under cases/> <case name> [y0 y1]"""
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402
import region_eval as RE                                             # noqa: E402
import region_refs as RR                                             # noqa: E402

B = VB.B
J = '/share/home/dq076/data/methane/Johnson2022_lake/Lake_Area.nc'


def main():
    ver, name = sys.argv[1], sys.argv[2]
    y0, y1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2010, 2019)
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, _, A_lake = VB.load_areas(case, years[0])
    e, parts = [], {k: [] for k in ('prod_tot', 'oxid_tot', 'surf_ebul', 'surf_diff')}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        e.append(((tr['f_methane_surf_flux_lake'].fillna(0) * A_land * w).sum('time') * B.M_CH4).values)
        for k in parts:                        # set C: per lake area, weighted by the lake area
            parts[k].append(((tr[f'f_methane_{k}_lake'].fillna(0) * A_lake.fillna(0) * w).sum('time')
                             * B.M_CH4).values)
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    E = np.mean(e, 0)
    PT = {k: np.mean(v, 0) for k, v in parts.items()}
    AL = A_lake.fillna(0).values / 1e12                                    # Mkm2
    names, W = RE.region_weights(lat, lon)
    with nc.Dataset(J) as h:
        jlat, jlon = h['lat'][:], h['lon'][:]
        ja = np.nan_to_num(np.ma.filled(h['LakeArea_Total'][:].astype(float), 0.0)) / 1e12
    jn, jr = RR.mask_on(jlat, jlon)
    rf = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                                  'region_refs.csv')).set_index('region')
    print(f'# {ver}/{name}, {years[0]}-{years[-1]}; lake area Mkm2, CH4 Tg/yr, per area g CH4 m-2 yr-1')
    print(f'{"region":24s} {"mod A":>6s} {"J A":>6s} {"mod E":>6s} {"J E":>6s} {"mod E/A":>7s} {"J E/A":>6s} | '
          f'per model lake area: {"P":>6s} {"O":>6s} {"ebul":>6s} {"diff":>6s}')
    tot = np.zeros(4)
    for r, n in enumerate(names):
        a, ea = np.nansum(W[r] * AL), np.nansum(W[r] * E)
        b, eb = ja[jr == r].sum(), rf.loc[n, 'lake_johnson']
        tot += (a, b, ea, eb)
        pp = {k: np.nansum(W[r] * v) / max(a * 1e12, 1e-9) * 1e12 for k, v in PT.items()}
        print(f'{n:24s} {a:6.3f} {b:6.3f} {ea:6.1f} {eb:6.1f} {ea / max(a, 1e-9):7.1f} {eb / max(b, 1e-9):6.1f} | '
              f'{"":20s} {pp["prod_tot"]:6.1f} {pp["oxid_tot"]:6.1f} {pp["surf_ebul"]:6.1f} {pp["surf_diff"]:6.1f}')
    print(f'{"total":24s} {tot[0]:6.3f} {tot[1]:6.3f} {tot[2]:6.1f} {tot[3]:6.1f} '
          f'{tot[2] / tot[0]:7.1f} {tot[3] / tot[1]:6.1f}')


if __name__ == '__main__':
    main()
