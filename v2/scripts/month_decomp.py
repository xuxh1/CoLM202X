#!/usr/bin/env python
"""Where the seasonal amplitude of tower CH4 goes missing: monthly climatology
per tower of the tower CH4 emission E_obs and ecosystem respiration Reco
(gap-filled, months with at least 80% of the steps), and of the model
heterotrophic respiration HR, CH4 production P, oxidation O and emission E,
the ratios E/(Reco/2) (tower) and E/HR, P/HR, O/P (model), and the 10 cm soil
temperature of tower and model. Fluxes in mg C m-2 d-1, ratios in %. The
last line gives the apparent Q10 (fit of ln flux on the 13 cm soil temperature
over the climatological months above 1 degC): tower E against the tower TS_1
(model temperature where the tower has none), model E, HR, P and O against
the model temperature.
Usage: month_decomp.py <case under cases/> <tag> site [site ...]"""
import glob
import sys
from collections import defaultdict

import netCDF4 as nc
import numpy as np

C = '/share/home/dq076/mode/Methane/cases/'
OBS = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Observation/'
DAY = 86400.


def obs(s):
    acc = {k: defaultdict(list) for k in ('E', 'R', 'T')}
    for p in glob.glob(f'{OBS}{s}_*_Flux.nc'):
        with nc.Dataset(p) as d:
            t = nc.num2date(d['time'][:], d['time'].units)
            mo = np.array([x.month for x in t])
            for k, v, fac in (('E', 'FCH4_f', 1e-9 * 12e3 * DAY), ('R', 'Resp_DT', 1e-6 * 12e3 * DAY), ('T', 'TS_1', 1.)):
                if v not in d.variables:
                    continue
                a = np.ma.filled(d[v][:].astype(float), np.nan).ravel() * fac
                if k == 'T' and np.nanmean(a) > 200:
                    a = a - 273.15
                for m in range(1, 13):
                    x = a[mo == m]
                    if x.size and np.isfinite(x).mean() >= 0.8:
                        acc[k][m].append(np.nanmean(x))
    return {k: np.array([np.mean(acc[k][m]) if acc[k][m] else np.nan for m in range(1, 13)]) for k in acc}


def model(case, tag, s):
    base = f'{C}{case}/sites_{tag}/{s}/history/{s}'
    acc = {k: [[] for _ in range(12)] for k in ('HR', 'P', 'O', 'E', 'T')}
    for f in sorted(glob.glob(base + '_hist_20??.nc')):
        y = f[-7:-3]
        with nc.Dataset(f) as d, nc.Dataset(base + f'_hist_tracer_{y}.nc') as t:
            n = d['f_hr'].shape[0]
            hr = np.ma.filled(d['f_hr'][:].astype(float), np.nan).reshape(n, -1)[:, 0] * 1e3 * DAY
            ts = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan).reshape(n, -1, d['f_t_soisno'].shape[-1])[:, 0, :]
            lev = list(np.asarray(d['soilsnow'][:]).astype(int))
            t10 = ts[:, lev.index(4)] - 273.15                   # soil layer 4, 0.09-0.17 m
            # a tower's first year may start mid-year (ID-Pag: June); records
            # carry mid-month time stamps
            mon = [x.month - 1 for x in nc.num2date(d['time'][:], d['time'].units)]
            for k, v in (('P', 'f_methane_prod_tot'), ('O', 'f_methane_oxid_tot'), ('E', 'f_methane_surf_flux_tot_active')):
                x = np.ma.filled(t[v][:].astype(float), np.nan).reshape(n, -1)[:, 0] * 12e3 * DAY
                for i, m in enumerate(mon):
                    acc[k][m].append(x[i])
            for i, m in enumerate(mon):
                acc['HR'][m].append(hr[i])
                acc['T'][m].append(t10[i])
    return {k: np.array([np.nanmean(v[m]) if v[m] else np.nan for m in range(12)]) for k, v in acc.items()}


def q10(f, t):
    ok = np.isfinite(f) & np.isfinite(t) & (t > 1.) & (f > 0)
    if ok.sum() < 4:
        return np.nan
    return float(np.exp(10. * np.polyfit(t[ok], np.log(f[ok]), 1)[0]))


def row(name, a, fmt='{:6.0f}'):
    return f'  {name:14s}' + ' '.join(fmt.format(x) if np.isfinite(x) else '     .' for x in a)


def main():
    case, tag = sys.argv[1], sys.argv[2]
    for s in sys.argv[3:]:
        o, m = obs(s), model(case, tag, s)
        print(f'{s}: month      ' + ' '.join(f'{k:6d}' for k in range(1, 13)))
        print(row('E obs', o['E']))
        print(row('E model', m['E']))
        print(row('Reco/2 obs', o['R'] / 2))
        print(row('HR model', m['HR']))
        print(row('E/(Reco/2) %', 100 * o['E'] / (o['R'] / 2), '{:6.1f}'))
        print(row('E/HR model %', 100 * m['E'] / m['HR'], '{:6.1f}'))
        print(row('P/HR model %', 100 * m['P'] / m['HR'], '{:6.1f}'))
        print(row('O/P model %', 100 * m['O'] / m['P'], '{:6.0f}'))
        print(row('T obs TS_1', o['T'], '{:6.1f}'))
        print(row('T model 13cm', m['T'], '{:6.1f}'))
        to = np.where(np.isfinite(o['T']), o['T'], m['T'])
        print('  apparent Q10  E obs {:.2f} | model E {:.2f} HR {:.2f} P {:.2f} O {:.2f}'.format(
            q10(o['E'], to), *(q10(m[k], m['T']) for k in ('E', 'HR', 'P', 'O'))))


if __name__ == '__main__':
    main()
