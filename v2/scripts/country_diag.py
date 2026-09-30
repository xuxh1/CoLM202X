#!/usr/bin/env python
"""Area or intensity: per country, CH4 (Tg/yr) and area (Mkm2, per-cell annual
maximum) of the permanent-wetland tile, the floodplain (as flood_regions.py)
and the lakes of a V2 global case, against WAD2M v2.0 annual-max area
(2010-2019), the GCP wetland 2025 22-run median (country_refs.csv) and
Johnson et al. (2022) lake area and CH4. Intensities in g CH4 m-2 yr-1.
Usage: country_diag.py <version dir under cases/> <case name> y0 y1 ISO [ISO ...]"""
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402
import country_eval as CE                                            # noqa: E402
import region_refs as RR                                             # noqa: E402

B = VB.B
V2 = CE.V2
J = '/share/home/dq076/data/methane/Johnson2022_lake'


def to_country(field, flat, flon):
    m = np.load(f'{V2}/data/country_mask_0p1.npz')
    g, isos = m['grid'], list(m['isos'])
    ii = np.clip(np.round((89.95 - np.asarray(flat)) / 0.1).astype(int), 0, 1799)
    fl = np.where(np.asarray(flon) > 180, np.asarray(flon) - 360, np.asarray(flon))
    jj = np.clip(np.round((fl + 179.95) / 0.1).astype(int), 0, 3599)
    out = np.bincount(g[np.ix_(ii, jj)].ravel(), weights=np.nan_to_num(field).ravel(),
                      minlength=len(isos) + 1)
    return dict(zip(isos, out[1:]))


def main():
    ver, name, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    want = sys.argv[5:]
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    A_land, _, A_lake = VB.load_areas(case, years[0])
    acc = {k: [] for k in ('et', 'ef', 'el', 'at', 'af')}
    for y in years:
        tr = xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc')
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')

        def tg(v):
            return ((tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4).values

        fsoil = tr['f_methane_soil_finundated'].fillna(0)
        dry = (fsoil.max('time') < 1e-4).values
        acc['et'].append(tg('f_methane_surf_flux_wetland'))
        acc['ef'].append(np.where(dry, 0, tg('f_methane_surf_flux_soil')))
        acc['el'].append(tg('f_methane_surf_flux_lake'))
        acc['at'].append((tr['f_methane_area_wetland'].fillna(0).max('time') * A_land / 1e12).values)
        acc['af'].append(((fsoil * tr['f_methane_area_soil'].fillna(0)).max('time') * A_land / 1e12).values)
        lat, lon = tr['lat'].values, tr['lon'].values
        tr.close()
    m = {k: np.mean(v, 0) for k, v in acc.items()}
    m['al'] = np.nan_to_num(A_lake.values) / 1e12
    isos, W = CE.country_weights(lat, lon)
    mod = {k: {iso: np.nansum(W[c + 1] * np.nan_to_num(v)) for c, iso in enumerate(isos)} for k, v in m.items()}
    la, lo, wa = RR.wad2m()
    wad = to_country(wa, la, lo)
    with nc.Dataset(f'{J}/Lake_Area.nc') as h:
        jla, jlo = h['lat'][:], h['lon'][:]
        a = np.nan_to_num(np.ma.filled(h['LakeArea_Total'][:].astype(np.float64), 0.0)) / 1e12
    jarea = to_country(a, jla, jlo)
    rf = pd.read_csv(f'{V2}/data/country_refs.csv').set_index('iso')
    print(f'# {ver}/{name} {years[0]}-{years[-1]}; E Tg/yr, A Mkm2, E/A g m-2 yr-1')
    print(f'{"iso":4s} {"Etile":>6s} {"Efl":>6s} {"Ewet":>6s} {"GCPmed":>6s} | {"Atile":>6s} {"Afl":>6s} {"WAD2M":>6s} | '
          f'{"tile/A":>6s} {"fl/A":>6s} {"GCP/WAD":>7s} | {"Elake":>6s} {"EJohn":>6s} {"Alake":>6s} {"AJohn":>8s} {"lk/A":>5s} {"J/A":>5s}')
    for iso in want:
        x = {k: mod[k][iso] for k in mod}
        r = rf.loc[iso]
        print(f'{iso:4s} {x["et"]:6.2f} {x["ef"]:6.2f} {x["et"] + x["ef"]:6.2f} {r.wet_med:6.2f} | '
              f'{x["at"]:6.3f} {x["af"]:6.3f} {wad[iso]:6.3f} | '
              f'{x["et"] / max(x["at"], 1e-9):6.1f} {x["ef"] / max(x["af"], 1e-9):6.1f} {r.wet_med / max(wad[iso], 1e-9):7.1f} | '
              f'{x["el"]:6.2f} {r.lake_johnson:6.2f} {x["al"]:6.3f} {jarea[iso]:8.3f} {x["el"] / max(x["al"], 1e-9):5.1f} {r.lake_johnson / max(jarea[iso], 1e-9):5.1f}')


if __name__ == '__main__':
    main()
