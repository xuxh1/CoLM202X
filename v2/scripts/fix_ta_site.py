#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Replace a tower's reanalysis air temperature in its Met forcing with the
measured air temperature of a nearby tower (user 2026-09-30, Q-66).

FLUXNET-CH4 has no measured TA at US-A10: its TA_F is ERA-Interim throughout and
runs 1.8-5.3 C warmer than US-NGB 5 km away in every summer month, while US-Beo
and US-Bes agree with US-NGB within 0.6 C (v2/results/diag_tundra_high_260930.txt).
The target's Tair becomes the donor's TA_F + 273.15 at the same half-hour stamp.
Its Qair is recomputed with the build rule of data/FLUXNET-CH4/build_colm_dataset.py
(es over water, Magnus; vapour pressure floored at 1 % RH) from the donor's RH_F at
the new temperature, capped at saturation, as the US-Beo fix took US-Bes humidity
(README section 6): the US-A10 RH record reads 27-60 % in every month against
79-90 % at US-NGB (v2/_scratch/forcing_ta/check_rh.py); half-hours without a donor
RH_F keep their vapour pressure. Psurf and the other variables are left alone. The Met file
is rewritten in place (the data set keeps one copy, data/FLUXNET-CH4/README.md
section 6); its md5 before the change goes to data/_tmp/分析_260930_塔气温与作物/.
Usage: fix_ta_site.py <target> <donor> [--apply]   (run on a compute node)"""
import datetime as dt
import glob
import hashlib
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

D = '/share/home/dq076/mode/Methane/data'
RAW = f'{D}/sites/FLUXNET-CH4/raw'
MET = f'{D}/FLUXNET-CH4/Forcing'
LOG = f'{D}/_tmp/分析_260930_塔气温与作物'
RH_FLOOR = 0.01


def raw(site, cols):
    f = glob.glob(f'{RAW}/FLX_{site}_FLUXNET-CH4_*/FLX_{site}_FLUXNET-CH4_HH_*.csv')
    assert len(f) == 1, f
    x = pd.read_csv(f[0], usecols=['TIMESTAMP_START'] + cols, na_values=[-9999])
    x.index = pd.to_datetime(x['TIMESTAMP_START'].astype(str), format='%Y%m%d%H%M')
    return x[cols]


def es_hpa(tc):
    return 6.1078 * 10 ** (7.5 * tc / (tc + 237.3))


def main():
    tgt, don = sys.argv[1], sys.argv[2]
    apply = '--apply' in sys.argv[3:]
    p = glob.glob(f'{MET}/{tgt}_*_Met.nc')
    assert len(p) == 1, p
    p = p[0]
    t = raw(tgt, ['TA_F'])
    dn = raw(don, ['TA_F', 'RH_F'])
    with nc.Dataset(p) as d:
        ta_old = d['Tair'][:, 0, 0].filled(np.nan)
        q_old = d['Qair'][:, 0, 0].filled(np.nan)
        ps = d['Psurf'][:, 0, 0].filled(np.nan)
    assert len(ta_old) == len(t), (len(ta_old), len(t))
    dev = np.nanmax(np.abs(ta_old - (t['TA_F'].values + 273.15)))
    print(f'{tgt}: {os.path.basename(p)}, {len(t)} half-hours; Met Tair vs raw TA_F max |diff| {dev:.4f} K')
    assert dev < 0.01, 'Met rows do not line up with the raw half-hours'
    tnew_c = dn['TA_F'].reindex(t.index).values
    have = np.isfinite(tnew_c)
    print(f'donor {don}: TA_F at {have.mean():.4f} of the target half-hours')
    ta_new = np.where(have, tnew_c + 273.15, ta_old)
    phpa = ps / 100.0
    ea_old = q_old * phpa / (0.622 + 0.378 * q_old)
    tc = ta_new - 273.15
    es = es_hpa(tc)
    rh = dn['RH_F'].reindex(t.index).values / 100.0
    ea = np.where(have & np.isfinite(rh), rh * es, ea_old)
    ea = np.clip(ea, RH_FLOOR * es, es)
    q_new = 0.622 * ea / (phpa - 0.378 * ea)
    m = t.index.month
    tab = pd.DataFrame({'dT': ta_new - ta_old, 'Tnew': ta_new - 273.15, 'Told': ta_old - 273.15,
                        'q_old': q_old * 1e3, 'q_new': q_new * 1e3}).groupby(m).mean()
    print('month  Told  Tnew    dT  q_old  q_new (g/kg)')
    for mm, r in tab.iterrows():
        print(f'{mm:5d} {r.Told:5.1f} {r.Tnew:5.1f} {r.dT:+5.2f}  {r.q_old:5.3f}  {r.q_new:5.3f}')
    if not apply:
        print('dry run: nothing written (add --apply)')
        return
    os.makedirs(LOG, exist_ok=True)
    with open(f'{LOG}/md5_before_forcing.txt', 'a') as f:
        f.write(f'{hashlib.md5(open(p, "rb").read()).hexdigest()}  {p}\n')
    with nc.Dataset(p, 'a') as d:
        d['Tair'][:, 0, 0] = ta_new
        d['Qair'][:, 0, 0] = q_new
        stamp = dt.datetime.now().strftime('%Y-%m-%d %H:%M')
        d.history = (getattr(d, 'history', '') + f'\n{stamp} V2 Q-66 fix_ta_site.py: Tair = {don} TA_F '
                     f'(measured, 5 km) at the same half-hour where available ({have.mean():.4f}), '
                     f'Qair recomputed from {don} RH_F at the new Tair (build rule, capped at saturation)').strip()
    print(f'written; md5 before in {LOG}/md5_before_forcing.txt')


if __name__ == '__main__':
    main()
