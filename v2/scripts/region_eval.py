#!/usr/bin/env python
"""Regional check of a V2 global case against the Global Methane Budget, by
the 18 GCP land regions (GCP2023_regions_mask.nc, the mask Saunois et al. 2025
use; USA, Canada, Brazil, Russia and China are single countries).

Model, per region, mean over the years given:
  wetland  = wetland tile + floodplain (soil tile in cells whose soil tile
             floods in the year), as global_eval.py and v2/目标与判据.md
  rice, lake = the source-split surface fluxes, land-area weighted
  area     = inundated area, sum of per-cell annual maxima (Mkm2)
Reference, 2010-2019, from Global_methane_budget_2023_2000_2020_v1.xlsx:
  wetland BU = 22 prognostic model runs, TD = 26 inversions: median [min-max]
  rice       = 7 inventories (CEDS, EDGARv6/v7, FAO, GAINS x2, USEPA): median [min-max]
  freshwater = lakes, ponds, reservoirs and rivers (Johnson or Lauerwald
               lakes + Rocher-Ros rivers): min of the two mins, mean of the
               two means, max of the two maxes; the model has lakes only.
Reference from v2/data/region_refs.csv (region_refs.py):
  lakeJ      = Johnson et al. (2022) lakes only, 41.6 Tg/yr globally
  WAD2M      = WAD2M v2.0 area, 2010-2019 mean of per-pixel annual maxima
A 2-degree cell is split between regions by the share of its four 1-degree
mask cells in each land region; a cell whose four mask cells are all ocean
goes to the nearest land mask cell.
Usage: region_eval.py <version dir under cases/> <case name> [y0 y1]"""
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
G = '/share/home/dq076/data/methane/GCP/GCP-CH4-2024'
XL = f'{G}/Global_methane_budget_2023_2000_2020_v1.xlsx'
NREG = 18                                   # region 18 is 'Oceans'


def gcp_table(sheet, category):
    """Rows of one category of a 2010-2019 budget sheet, columns = region names."""
    df = pd.read_excel(XL, sheet_name=sheet, header=None)
    h = df.index[df[0].astype(str).str.strip() == 'category'][0]
    cols = [str(c) for c in df.iloc[h].values]
    body = df.iloc[h + 1:].copy()
    body.columns = cols
    body['category'] = body['category'].ffill()
    src = cols[1]
    sub = body[body['category'].astype(str).str.strip() == category].set_index(src)
    return sub.drop(columns=['category']).apply(pd.to_numeric, errors='coerce')


def region_weights(lat, lon):
    """(nreg, nlat, nlon) share of each 2-degree cell in each land region."""
    m = xr.open_dataset(f'{G}/GCP2023_regions_mask.nc')
    names = [str(s).strip().replace(' ', '_') for s in m['regions_name'].values]
    reg = m['regions'].values.astype(int)
    mlat, mlon = m['lat'].values, m['lon'].values
    mlon = np.where(mlon > 180, mlon - 360, mlon)
    land = np.argwhere(reg != NREG)
    W = np.zeros((NREG, lat.size, lon.size))
    dlat = abs(lat[1] - lat[0]) / 2
    dlon = abs(lon[1] - lon[0]) / 2
    for i, la in enumerate(lat):
        ii = np.where(np.abs(mlat - la) < dlat)[0]
        for j, lo in enumerate(lon):
            lo = lo - 360 if lo > 180 else lo
            jj = np.where(np.abs(((mlon - lo) + 180) % 360 - 180) < dlon)[0]
            r = reg[np.ix_(ii, jj)].ravel()
            r = r[r != NREG]
            if r.size == 0:
                # nearest land mask cell (great-circle distance is overkill here)
                d = (mlat[land[:, 0]] - la) ** 2 + ((((mlon[land[:, 1]] - lo) + 180) % 360 - 180)
                                                    * np.cos(np.radians(la))) ** 2
                r = np.array([reg[tuple(land[np.argmin(d)])]])
            for k in r:
                W[k, i, j] += 1.0 / r.size
    return names[:NREG], W


def main():
    ver, name = sys.argv[1], sys.argv[2]
    y0, y1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2010, 2019)
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, A_active, A_lake = VB.load_areas(case, years[0])
    acc = {k: [] for k in ('wet', 'rice', 'lake', 'area')}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')

        def tg(v):
            return ((tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4).values

        es = tg('f_methane_surf_flux_soil')
        dry = (tr['f_methane_soil_finundated'].fillna(0).max('time') < 1e-4).values
        acc['wet'].append(tg('f_methane_surf_flux_wetland') + np.where(dry, 0, es))
        acc['rice'].append(tg('f_methane_surf_flux_rice'))
        acc['lake'].append(tg('f_methane_surf_flux_lake'))
        fin = tr['f_methane_finundated'].fillna(0)
        acc['area'].append(((fin * A_active).max('time') / 1e12).fillna(0).values)
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    mean = {k: np.nanmean(np.array(v), 0) for k, v in acc.items()}
    names, W = region_weights(lat, lon)
    reg = {k: np.array([np.nansum(W[r] * mean[k]) for r in range(NREG)]) for k in mean}

    bu = gcp_table('2010-2019 BU Budget', 'Wetlands')
    td = gcp_table('2010-2019 TD Budget', 'Wetlands')
    ri = gcp_table('2010-2019 BU Budget', 'Rice cultivation')
    fw = gcp_table('2010-2019 BU Budget', 'Freshwaters')
    rf = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                                  'region_refs.csv')).set_index('region')

    def stat(t, c):
        v = t[c].dropna().values
        return np.median(v), v.min(), v.max()

    def flag(x, lo, hi):
        return 'in ' if lo <= x <= hi else ('LOW' if x < lo else 'HI ')

    L = [f'# {ver}/{name}, {years[0]}-{years[-1]} ({len(years)} years); Tg CH4/yr; area Mkm2',
         '# wetland = wetland tile + floodplain; BU 22 model runs, TD 26 inversions, rice 7 inventories: median [min-max]',
         f'{"region":24s} {"wetland":>7s} {"BU med [min-max]":>20s} bu  {"TD med [min-max]":>20s} td '
         f'| {"rice":>5s} {"inv med [min-max]":>18s} ri | {"lake":>5s} {"freshw mean [min-max]":>21s} {"lakeJ":>5s} | {"area":>5s} {"WAD2M":>5s}']
    tot = {k: 0.0 for k in reg}
    for r, n in enumerate(names):
        b, t, q = stat(bu, n), stat(td, n), stat(ri, n)
        fmean = np.nanmean([fw.loc[i, n] for i in fw.index if i.endswith('-mean')])
        fmin = np.nanmin([fw.loc[i, n] for i in fw.index if i.endswith('-min')])
        fmax = np.nanmax([fw.loc[i, n] for i in fw.index if i.endswith('-max')])
        x = {k: reg[k][r] for k in reg}
        for k in tot:
            tot[k] += x[k]
        L.append(f'{n:24s} {x["wet"]:7.1f} {b[0]:6.1f} [{b[1]:5.1f}-{b[2]:5.1f}] {flag(x["wet"], b[1], b[2])} '
                 f'{t[0]:6.1f} [{t[1]:5.1f}-{t[2]:5.1f}] {flag(x["wet"], t[1], t[2])} '
                 f'| {x["rice"]:5.1f} {q[0]:5.1f} [{q[1]:4.1f}-{q[2]:5.1f}] {flag(x["rice"], q[1], q[2])} '
                 f'| {x["lake"]:5.1f} {fmean:5.1f} [{fmin:5.1f}-{fmax:5.1f}] {rf.loc[n, "lake_johnson"]:5.1f} '
                 f'| {x["area"]:5.2f} {rf.loc[n, "wad2m_area"]:5.2f}')
    b, t, q = stat(bu, 'Global'), stat(td, 'Global'), stat(ri, 'Global')
    L.append(f'{"sum of land regions":24s} {tot["wet"]:7.1f} {b[0]:6.1f} [{b[1]:5.1f}-{b[2]:5.1f}]     '
             f'{t[0]:6.1f} [{t[1]:5.1f}-{t[2]:5.1f}]     | {tot["rice"]:5.1f} {q[0]:5.1f} [{q[1]:4.1f}-{q[2]:5.1f}]    '
             f'| {tot["lake"]:5.1f}                       {rf["lake_johnson"].sum():5.1f} '
             f'| {tot["area"]:5.2f} {rf["wad2m_area"].sum():5.2f}')
    nb = sum(stat(bu, n)[1] <= reg['wet'][r] <= stat(bu, n)[2] for r, n in enumerate(names))
    nt = sum(stat(td, n)[1] <= reg['wet'][r] <= stat(td, n)[2] for r, n in enumerate(names))
    nr = sum(stat(ri, n)[1] <= reg['rice'][r] <= stat(ri, n)[2] for r, n in enumerate(names))
    L.append(f'# regions inside range: wetland BU {nb}/18, TD {nt}/18; rice {nr}/18')
    # log distance to the BU median, a single number to compare runs by
    lb = [abs(np.log(max(reg['wet'][r], 0.05) / max(stat(bu, n)[0], 0.05))) for r, n in enumerate(names)]
    L.append(f'# median |ln(model / BU median)| over regions: {np.median(lb):.2f}')
    ll = [abs(np.log(max(reg['lake'][r], 0.05) / rf.loc[n, 'lake_johnson'])) for r, n in enumerate(names)]
    la = [abs(np.log(max(reg['area'][r], 0.005) / rf.loc[n, 'wad2m_area'])) for r, n in enumerate(names)]
    L.append(f'# median |ln(model / Johnson)| lakes: {np.median(ll):.2f}; |ln(model / WAD2M)| area: {np.median(la):.2f}')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
