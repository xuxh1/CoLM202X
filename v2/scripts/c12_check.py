#!/usr/bin/env python
"""C-12 check (ledger 2.11): does heterotrophic respiration of the wetland patch
now follow productivity? Per tower, over months with observed ecosystem
respiration (monthly means need >= 480 valid half-hours): annual model HR under
V-2 (v2_c11_sp30, fixed substrate) and V-5 (v2_c12_sp30), the plant litter
input of V-5 (0.5 x biophysical photosynthesis), observed ER and GPP; cross-tower
Spearman of HR with observed ER; plus the CH4 scores of both tags."""
import csv
import glob

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
V2 = f'{R}/CoLM202X-paper-v2/v2'
OBS = f'{R}/data/FLUXNET-CH4/Observation'
RUNS = {'V2': ('paper_v2/v260924b/sp_v2_c11', 'v2_c11_sp30'), 'V5': ('paper_v2/v260924d/sp_v2_c12', 'v2_c12_sp30')}
SITES = [(r['ID'], int(r['SITE_wetland_class'])) for r in csv.DictReader(open(f'{R}/scripts/sites/LIST_sites_wtd23.csv'))]
CLS = {1: 'permafrost', 3: 'bog', 4: 'fen', 5: 'marsh', 6: 'salt_marsh', 7: 'trop_swamp'}
SALT = {'DE-Hte', 'US-Srr', 'US-StJ'}


def obs_monthly(site):
    f = sorted(glob.glob(f'{OBS}/{site}_*_Flux.nc'))[0]
    with nc.Dataset(f) as d:
        t = nc.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False,
                        only_use_python_datetimes=True)
        g = lambda v: np.ma.filled(d[v][:], np.nan).ravel().astype(float) if v in d.variables else np.full(len(t), np.nan)
        df = pd.DataFrame({'er': g('Resp'), 'gpp': g('GPP')}, index=pd.DatetimeIndex(t))
    cnt, mean = df.resample('MS').count(), df.resample('MS').mean()
    mean[cnt < 480] = np.nan
    return mean * 12.011e-6 * 86400 * mean.index.days_in_month.values[:, None]   # g C m-2 per month


def model(case, tag, site):
    rows = []
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{site}/history/{site}_hist_2*.nc')):
        y = int(f[-7:-3])
        with nc.Dataset(f) as h:
            g = lambda v: np.ma.filled(h[v][:], np.nan).astype(float).ravel()
            hr = g('f_hr')
            if hr.size != 12:
                continue
            a = g('f_assim')
            soc = g('f_totsomc') + g('f_totlitc')
            for i in range(12):
                sec = pd.Timestamp(y, i + 1, 1).days_in_month * 86400.0
                rows.append({'t': pd.Timestamp(y, i + 1, 1), 'hr': hr[i] * sec,
                             'inp': a[i] * 12.011 * 0.5 * sec, 'soc': soc[i]})
    return pd.DataFrame(rows).set_index('t')


def main():
    sc = {k: pd.read_csv(f'{V2}/results/scores/score_series_{tag}.csv').set_index('site') for k, (_, tag) in RUNS.items()}
    out = []
    for site, cls in SITES:
        o = obs_monthly(site)
        r = {'site': site, 'cls': CLS[cls]}
        for k, (case, tag) in RUNS.items():
            m = model(case, tag, site)
            j = m.join(o, how='inner').dropna(subset=['er'])
            r[f'hr_{k}'] = m.hr.sum() / (len(m) / 12)
            if k == 'V5':
                r['input_V5'] = m.inp.sum() / (len(m) / 12)
                r['soc_V5_kgm2'] = m.soc.mean() / 1000
            if len(j) >= 6:
                r[f'hr_over_obsER_{k}'] = j.hr.sum() / j.er.sum()
                r['obsER'] = j.er.mean() * 12
                r['obsGPP'] = j.gpp.mean() * 12 if j.gpp.notna().sum() >= 6 else np.nan
            for v in ('KGEln', 'beta', 'r'):
                r[f'{v}_{k}'] = sc[k].loc[site, v]
        out.append(r)
    t = pd.DataFrame(out)
    x = t[~t.site.isin(SALT)].dropna(subset=['obsER'])
    L = ['# C-12 check: annual HR, plant input (0.5 x assim) in g C m-2 yr-1; HR / observed ER over matched months',
         f'# {len(x)} non-saline towers with observed ER']
    for k in RUNS:
        q = np.log10(x[f'hr_over_obsER_{k}'])
        L.append(f'{k}: median HR/obsER {x[f"hr_over_obsER_{k}"].median():.2f}  log10 sd {q.std():.2f}  '
                 f'spearman HR vs obsER {x[f"hr_{k}"].corr(x.obsER, method="spearman"):.2f}')
    L.append(f'V5 input vs obsGPP: median ratio {(x.input_V5 / x.obsGPP).median():.2f}')
    L.append('')
    cols = ['site', 'cls', 'obsER', 'obsGPP', 'hr_V2', 'hr_V5', 'input_V5', 'soc_V5_kgm2',
            'KGEln_V2', 'KGEln_V5', 'beta_V2', 'beta_V5', 'r_V2', 'r_V5']
    L.append(t[cols].round(2).to_string(index=False))
    open(f'{V2}/results/c12_check.log', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
