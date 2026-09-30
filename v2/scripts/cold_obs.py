"""Cold-season check against the tower's own records (Q-27): per calendar
month, the share of half-hours with a measured CH4 flux, the mean of the
measured half-hours beside the gap-filled mean, the gap-filled mean of each
year for September to January, the tower soil temperature at its first three
probes (depths from the FLUXNET-CH4 META table, which writes metres with a
'cm' suffix) beside the model soil temperature interpolated to the same depth,
and tower snow depth beside the model's. Model months are the same years as
the tower record.
Usage: cold_obs.py <case under cases/> <tag> <site> [site ...]"""
import csv
import glob
import re
import sys
from collections import defaultdict

import netCDF4 as nc
import numpy as np

M = '/share/home/dq076/mode/Methane'
OBS = f'{M}/data/FLUXNET-CH4/Observation'
META = f'{M}/data/sites/FLUXNET-CH4/raw/FLX_AA-Flx_CH4-META_20201112135337801132.csv'
TO_MG = 16.043e-6 * 86400                # nmol m-2 s-1 -> mg CH4 m-2 d-1
ZSOI = 0.025 * (np.exp(0.5 * (np.arange(1, 11) - 0.5)) - 1.)     # CoLM node depths, m


def probes():
    out = {}
    with open(META, encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            out[r['SITE_ID']] = {int(k): abs(float(v)) for k, v in
                                 re.findall(r'TS_(\d+)\s*=\s*(-?[\d.]+)\s*cm', r['SOIL_TEMP_PROBE_DEPTHS'] or '')}
    return out


def tower(sid):
    p = glob.glob(f'{OBS}/{sid}_*_Flux.nc')[0]
    with nc.Dataset(p) as d:
        t = nc.num2date(d['time'][:], d['time'].units)
        get = lambda v: np.ma.filled(d[v][:].astype(float), np.nan).ravel() if v in d.variables else None
        return t, {v: get(v) for v in ('FCH4', 'FCH4_f', 'TS_1', 'TS_2', 'TS_3', 'SnowDepth')}


def main():
    case, tag = sys.argv[1:3]
    dep = probes()
    for sid in sys.argv[3:]:
        t, v = tower(sid)
        ym = np.array([(x.year, x.month) for x in t])
        yrs = sorted({y for y, m in ym[np.isfinite(v['FCH4_f'])]})
        cov, meas, fill, per = {}, {}, {}, defaultdict(dict)
        ts = {k: {} for k in ('TS_1', 'TS_2', 'TS_3', 'SnowDepth')}
        for m in range(1, 13):
            sel = ym[:, 1] == m
            f = v['FCH4'][sel]
            cov[m] = np.isfinite(f).mean() if sel.any() else np.nan
            meas[m] = np.nanmean(f) * TO_MG if np.isfinite(f).any() else np.nan
            fill[m] = np.nanmean(v['FCH4_f'][sel]) * TO_MG if sel.any() else np.nan
            for k in ts:
                if v[k] is not None and np.isfinite(v[k][sel]).any():
                    ts[k][m] = np.nanmean(v[k][sel])
            for y in yrs:
                s2 = sel & (ym[:, 0] == y)
                if s2.any() and np.isfinite(v['FCH4_f'][s2]).any():
                    per[m][y] = np.nanmean(v['FCH4_f'][s2]) * TO_MG
        # model soil temperature at the probe depths and snow depth, same years
        mt = defaultdict(lambda: defaultdict(list))
        for y in yrs:
            f = f'{M}/cases/{case}/sites_{tag}/{sid}/history/{sid}_hist_{y}.nc'
            try:
                d = nc.Dataset(f)
            except OSError:
                continue
            with d:
                T = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)[:, 0, 5:15] - 273.15
                sd = np.ma.filled(d['f_snowdp'][:].astype(float), np.nan).reshape(T.shape[0], -1)[:, 0]
                for m in range(min(12, T.shape[0])):
                    for k, z in dep.get(sid, {}).items():
                        if k <= 3:
                            mt[f'TS_{k}'][m + 1].append(np.interp(z, ZSOI, T[m]))
                    mt['SnowDepth'][m + 1].append(sd[m])
        row = lambda xs, fmt='6.1f': ' '.join(f'{x:{fmt}}' if np.isfinite(x) else '     .' for x in xs)
        print(f'== {sid}  years {yrs[0]}-{yrs[-1]}  probes {dep.get(sid)}')
        print('  month               ' + ' '.join(f'{m:6d}' for m in range(1, 13)))
        print('  measured share      ' + row([cov[m] for m in range(1, 13)], '6.2f'))
        print('  FCH4 measured mean  ' + row([meas[m] for m in range(1, 13)]))
        print('  FCH4 gap-filled     ' + row([fill[m] for m in range(1, 13)]))
        for y in yrs:
            print(f'    {y} gap-filled    ' + row([per[m].get(y, np.nan) for m in range(1, 13)]))
        for k in ('TS_1', 'TS_2', 'TS_3', 'SnowDepth'):
            if ts[k]:
                fmt = '6.2f' if k == 'SnowDepth' else '6.1f'
                print(f'  tower {k:12s}  ' + row([ts[k].get(m, np.nan) for m in range(1, 13)], fmt))
                if mt[k]:
                    print(f'  model {k:12s}  ' + row([np.mean(mt[k][m]) if mt[k][m] else np.nan for m in range(1, 13)], fmt))


if __name__ == '__main__':
    main()
