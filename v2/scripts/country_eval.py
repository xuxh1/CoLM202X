#!/usr/bin/env python
"""Country check of a V2 global case (Tg CH4/yr, mean over the years given):
model wetland (wetland tile + floodplain, as region_eval.py), rice and lakes
per country against v2/data/country_refs.csv (country_refs.py):
  wetland  GCP-CH4 wetland synthesis 2025, 22 prognostic runs, 2010-2019:
           median [min-max]
  rice     GRPI (2022) and the GCP 2024 inversion prior (2010-2019); the flag
           is against [min-max] of the two
  lake     Johnson et al. (2022) lakes and the GCP 2024 freshwater prior
Countries listed: reference wetland median, rice or lake of at least 0.5 Tg/yr.
A 2-degree model cell is split between countries by the share of its land
0.1-degree cells (country_mask_0p1.npz) in each country.
Usage: country_eval.py <version dir under cases/> <case name> [y0 y1]"""
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def country_weights(lat, lon):
    """(ncountry+1, nlat, nlon) share of the land of each model cell per country;
    index 0 (no country) is kept only where a cell has no land 0.1-degree cell."""
    m = np.load(f'{V2}/data/country_mask_0p1.npz')
    g, glat, glon, isos = m['grid'], m['lat'], m['lon'], list(m['isos'])
    n = len(isos)
    W = np.zeros((n + 1, lat.size, lon.size))
    dlat, dlon = abs(lat[1] - lat[0]) / 2, abs(lon[1] - lon[0]) / 2
    for i, la in enumerate(lat):
        ii = np.where(np.abs(glat - la) < dlat)[0]
        wl = np.cos(np.radians(glat[ii]))
        for j, lo in enumerate(lon):
            lo = lo - 360 if lo > 180 else lo
            jj = np.where(np.abs(((glon - lo) + 180) % 360 - 180) < dlon)[0]
            c = g[np.ix_(ii, jj)]
            w = np.broadcast_to(wl[:, None], c.shape)
            land = c > 0
            if land.any():
                W[:, i, j] = np.bincount(c[land], weights=w[land], minlength=n + 1) / w[land].sum()
            else:
                W[0, i, j] = 1.0
    return isos, W


def main():
    ver, name = sys.argv[1], sys.argv[2]
    y0, y1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (2010, 2019)
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, _, _ = VB.load_areas(case, years[0])
    acc = {k: [] for k in ('wet', 'rice', 'lake')}
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
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    mean = {k: np.nanmean(np.array(v), 0) for k, v in acc.items()}
    isos, W = country_weights(lat, lon)
    mod = {k: np.array([np.nansum(W[c] * np.nan_to_num(mean[k])) for c in range(1, len(isos) + 1)])
           for k in mean}
    rf = pd.read_csv(f'{V2}/data/country_refs.csv').set_index('iso').loc[isos]
    for k in mod:
        rf[f'm_{k}'] = mod[k]
    keep = (rf.wet_med >= 0.5) | (rf[['rice_grpi', 'rice_prior']].max(1) >= 0.5) | \
           (rf[['lake_johnson', 'fresh_prior']].max(1) >= 0.5)
    t = rf[keep].sort_values('wet_med', ascending=False)

    def flag(x, lo, hi):
        return 'in ' if lo <= x <= hi else ('LOW' if x < lo else 'HI ')

    L = [f'# {ver}/{name}, {years[0]}-{years[-1]} ({len(years)} years); Tg CH4/yr; '
         f'{len(t)} countries with a reference wetland, rice or lake source of at least 0.5',
         '# wetland: GCP wetland 2025, 22 runs, median [min-max]; rice: GRPI / GCP prior; lake: Johnson 2022 / GCP freshwater prior',
         f'{"iso":4s} {"wetland":>7s} {"ref med [min-max]":>20s} fl | {"rice":>5s} {"GRPI":>5s} {"prior":>5s} fl | '
         f'{"lake":>5s} {"John":>5s} {"prior":>5s} fl']
    nw = nr = nl = cw = cr = cl = 0
    dw = []
    for iso, r in t.iterrows():
        fw = flag(r.m_wet, r.wet_min, r.wet_max) if r.wet_med >= 0.5 else '   '
        rlo, rhi = min(r.rice_grpi, r.rice_prior), max(r.rice_grpi, r.rice_prior)
        fr = flag(r.m_rice, rlo, rhi) if rhi >= 0.5 else '   '
        llo, lhi = min(r.lake_johnson, r.fresh_prior), max(r.lake_johnson, r.fresh_prior)
        fl = flag(r.m_lake, llo, lhi) if lhi >= 0.5 else '   '
        if r.wet_med >= 0.5:
            cw += 1
            nw += fw == 'in '
            dw.append(abs(np.log(max(r.m_wet, 0.05) / max(r.wet_med, 0.05))))
        if rhi >= 0.5:
            cr += 1
            nr += fr == 'in '
        if lhi >= 0.5:
            cl += 1
            nl += fl == 'in '
        L.append(f'{iso:4s} {r.m_wet:7.2f} {r.wet_med:6.2f} [{r.wet_min:5.2f}-{r.wet_max:5.2f}] {fw} | '
                 f'{r.m_rice:5.2f} {r.rice_grpi:5.2f} {r.rice_prior:5.2f} {fr} | '
                 f'{r.m_lake:5.2f} {r.lake_johnson:5.2f} {r.fresh_prior:5.2f} {fl}')
    L.append(f'# inside range: wetland {nw}/{cw} (median |ln(model/ref median)| {np.median(dw):.2f}); '
             f'rice {nr}/{cr}; lake {nl}/{cl}')
    L.append(f'# listed-country totals: model wetland {t.m_wet.sum():.1f} vs ref median {t.wet_med.sum():.1f}; '
             f'model rice {t.m_rice.sum():.1f}; model lake {t.m_lake.sum():.1f}; world model wetland {rf.m_wet.sum():.1f}')
    out = '\n'.join(L)
    print(out)
    os.makedirs(f'{V2}/results', exist_ok=True)
    with open(f'{V2}/results/country_{name}.txt', 'w') as f:
        f.write(out + '\n')


if __name__ == '__main__':
    main()
