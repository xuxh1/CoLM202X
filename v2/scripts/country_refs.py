#!/usr/bin/env python
"""Per-country reference CH4 (Tg CH4/yr, 2010-2019 unless the source is a
climatology), written once to v2/data/country_refs.csv for country_eval.py:
  wetland   GCP-CH4 wetland synthesis 2025 (Zhang et al. 2025,
            doi:10.5194/bg-22-305-2025), each prognostic run of
            netcdf_totflux.zip (kg CH4 m-2 s-1 per grid-cell area, monthly):
            2010-2019 mean per country, then median, min and max over runs
  rice      GRPI rice paddy CH4 (grpi_hemco.nc, 2022, 0.1 degree)
  lake      Johnson et al. (2022) lakes, climatology (as region_refs.py)
  rice_prior, fresh_prior   GCP 2024 inversion priors (Saunois et al. 2025,
            GCP_Prior_CH4_fluxes.nc, 1 degree): rice 2010-2019 mean,
            freshwaters climatology
GRPI has gaps (no rice in Colombia or Malaysia), hence the second rice source.
Countries: Natural Earth admin-0 polygons of ~/data/boundaries/world, ISO_A3,
rasterised at 0.1 degree (majority of the polygon covering each cell centre);
grids of the sources are mapped by nearest 0.1-degree cell.
Also writes v2/data/country_mask_0p1.npz (ISO index grid) for country_eval.py.
Usage: country_refs.py"""
import gzip
import os
import shutil
import tempfile

import geopandas as gpd
import netCDF4 as nc
import numpy as np
import pandas as pd
from rasterio import features
from rasterio.transform import from_origin

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
W = '/share/home/dq076/data/boundaries/world/world.shp'
GCPZ = '/share/home/dq076/data/methane/GCP/GCP-CH4-wetland-2025/netcdf_totflux.zip'
GRPI = '/share/home/dq076/data/methane/GRPI_rice/grpi_hemco.nc'
JOH = '/share/home/dq076/data/methane/Johnson2022_lake'
PRIOR = '/share/home/dq076/data/methane/GCP/GCP-CH4-2024/GCP_Prior_CH4_fluxes.nc'
R = 6371e3


def country_grid():
    w = gpd.read_file(W)
    w = w[w.ISO_A3_EH.str.len() == 3].copy()
    w['iso'] = w.ISO_A3_EH.where(w.ISO_A3_EH != '-99', w.ADM0_A3)
    isos = sorted(w.iso.unique())
    idx = {k: i + 1 for i, k in enumerate(isos)}
    shapes = [(g, idx[k]) for g, k in zip(w.geometry, w.iso)]
    tr = from_origin(-180, 90, 0.1, 0.1)
    grid = features.rasterize(shapes, out_shape=(1800, 3600), transform=tr, fill=0, dtype='int32')
    lat = 90 - 0.05 - 0.1 * np.arange(1800)
    lon = -180 + 0.05 + 0.1 * np.arange(3600)
    np.savez_compressed(f'{V2}/data/country_mask_0p1.npz', grid=grid, lat=lat, lon=lon, isos=np.array(isos))
    return grid, lat, lon, isos


def to_country(field, flat, flon, grid, lat, lon, n):
    """Sum a per-cell total (Tg) of a source grid into countries by the 0.1-degree
    country cell nearest to each source cell centre."""
    ii = np.clip(np.round((90 - 0.05 - np.asarray(flat)) / 0.1).astype(int), 0, 1799)
    fl = np.where(np.asarray(flon) > 180, np.asarray(flon) - 360, np.asarray(flon))
    jj = np.clip(np.round((fl + 180 - 0.05) / 0.1).astype(int), 0, 3599)
    cid = grid[np.ix_(ii, jj)]
    out = np.bincount(cid.ravel(), weights=np.nan_to_num(field).ravel(), minlength=n + 1)
    return out


def cell_area(lat, dlat, dlon):
    return (np.radians(dlat) * np.radians(dlon) * R ** 2 * np.cos(np.radians(lat)))[:, None]


def main():
    grid, lat, lon, isos = country_grid()
    n = len(isos)
    rows = {}
    tmp = tempfile.mkdtemp(dir='/tmp')
    import zipfile
    z = zipfile.ZipFile(GCPZ)
    runs = [m for m in z.namelist() if m.endswith('.nc.gz') and not m.startswith('__MACOSX')]
    for m in runs:
        name = os.path.basename(m).replace('_totflux', '').replace('_prognostic.nc.gz', '')
        gz = z.extract(m, tmp)
        f = gz[:-3]
        with gzip.open(gz, 'rb') as a, open(f, 'wb') as b:
            shutil.copyfileobj(a, b)
        os.remove(gz)
        with nc.Dataset(f) as d:
            la = d['lat'][:] if 'lat' in d.variables else d['latitude'][:]
            lo = d['lon'][:] if 'lon' in d.variables else d['longitude'][:]
            # every run is monthly from January 2000 (time units differ between
            # runs, several are non-CF), so the year is taken from the index
            yrs = 2000 + np.arange(d['time'].size) // 12
            sel = np.where((yrs >= 2010) & (yrs <= 2019))[0]
            x = np.nan_to_num(np.ma.filled(d['totflux'][sel].astype(np.float64), 0.0))
            x[x > 1e10] = 0
            mean = x.mean(0) * 86400 * 365.25                              # kg m-2 yr-1
        A = cell_area(np.asarray(la), abs(float(la[1] - la[0])), abs(float(lo[1] - lo[0])))
        tg = mean * A / 1e9
        rows[name] = to_country(tg, la, lo, grid, lat, lon, n)
        print(f'{name}: global {tg.sum():.1f} Tg/yr, years {yrs[sel][0]}-{yrs[sel][-1]}')
        os.remove(f)
    shutil.rmtree(tmp)
    wet = np.array(list(rows.values()))                                   # runs x countries
    # rice (GRPI 2022) and lakes (Johnson 2022)
    with nc.Dataset(GRPI) as d:
        la, lo = d['lat'][:], d['lon'][:]
        days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        A = cell_area(np.asarray(la), 0.1, 0.1)
        rice = np.zeros((la.size, lo.size))
        for k in range(12):
            rice += np.nan_to_num(np.ma.filled(d['emi_ch4'][k].astype(float), 0.0)) * A * days[k] * 86400 / 1e9
    rice_c = to_country(rice, la, lo, grid, lat, lon, n)
    with nc.Dataset(f'{JOH}/Lake_Area.nc') as h:
        la, lo = h['lat'][:], h['lon'][:]
    A = cell_area(np.asarray(la), 0.25, 0.25)
    lake = np.zeros((la.size, lo.size))
    for f, v in (('Lake_CH4_Diff_Ebul_Emiss.nc', 'DiffusionEbullition_TotalLakes'),
                 ('Lake_CH4_Ice_out_Emiss.nc', 'IceOut_TotalLakes'),
                 ('Lake_CH4_Fall_Turnover_Emiss.nc', 'FallTurnover_TotalLakes')):
        with nc.Dataset(f'{JOH}/{f}') as h:
            x = h[v]
            for t0 in range(0, x.shape[0], 30):
                a = np.nan_to_num(np.ma.filled(x[t0:t0 + 30].astype(np.float64), 0.0))
                lake += np.clip(a, 0, None).sum(0) * A / 1e12
    lake_c = to_country(lake, la, lo, grid, lat, lon, n)
    with nc.Dataset(PRIOR) as d:
        la, lo = d['lat'][:], d['lon'][:]
        A = cell_area(np.asarray(la), 1.0, 1.0)
        sec = 86400 * 365.25
        r = np.ma.filled(d['flux_ch4_rice'][120:240].astype(np.float64), 0.0).mean(0)
        fw = np.ma.filled(d['flux_ch4_freshwaters'][:].astype(np.float64), 0.0).mean(0)
    rp_c = to_country(np.nan_to_num(r) * A * sec / 1e9, la, lo, grid, lat, lon, n)
    fw_c = to_country(np.nan_to_num(fw) * A * sec / 1e9, la, lo, grid, lat, lon, n)
    df = pd.DataFrame({'iso': isos,
                       'wet_med': np.median(wet[:, 1:], 0), 'wet_min': wet[:, 1:].min(0),
                       'wet_max': wet[:, 1:].max(0), 'wet_n': wet.shape[0],
                       'rice_grpi': rice_c[1:], 'rice_prior': rp_c[1:],
                       'lake_johnson': lake_c[1:], 'fresh_prior': fw_c[1:]})
    df.to_csv(f'{V2}/data/country_refs.csv', index=False, float_format='%.4f')
    print(f'runs {wet.shape[0]}; global wetland median of run totals {np.median(wet.sum(1)):.1f}; '
          f'rice {rice_c.sum():.1f} (prior {rp_c.sum():.1f}); lakes {lake_c.sum():.1f} (freshwater prior {fw_c.sum():.1f})')
    print(df.sort_values('wet_med', ascending=False).head(20).to_string(index=False))


if __name__ == '__main__':
    main()
