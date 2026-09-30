#!/usr/bin/env python
"""Is a site tag's result settled? Compare it with the same configuration run
after a different spin-up length. Settled means more spin-up no longer changes
what is scored: the monthly CH4 flux and its four scores.

Per tower, against the reference tag:
  dE_y1, dE_last  relative change of the mean CH4 flux in the first and in the
                  last complete formal year (the first year replays the spin-up
                  year's forcing; the last shows whether a gap closes during
                  the formal run itself)
  dE_all          relative change of the mean flux over all formal months
  dSOM, dLIT      relative change of totsomc and totlitc in the first year
  dKGE, dbeta, dalpha, dr   score changes, when both score csv files exist
Only towers carrying a .complete marker in both tags are read.
Usage: eq_compare.py <case under cases/> <reference tag> <tag> [tag ...]"""
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
V2 = f'{R}/CoLM202X-paper-v2/v2'
CH4 = 'f_methane_surf_flux_tot_active'
SALT = {'DE-Hte', 'US-Srr', 'US-StJ'}


def annual(case, tag, site, var, kind):
    """{year: annual mean} over complete (12-record) history years."""
    pat = f'{R}/cases/{case}/sites_{tag}/{site}/history/{site}_hist_{kind}2*.nc'
    out = {}
    for f in sorted(glob.glob(pat)):
        with nc.Dataset(f) as h:
            if var not in h.variables:
                continue
            v = np.ma.filled(h[var][:], np.nan).reshape(h[var].shape[0], -1)[:, 0]
        if v.size == 12 and np.all(np.isfinite(v)) and np.all(np.abs(v) < 1e30):
            out[int(f[-7:-3])] = float(v.mean())
    return out


def rel(a, b):
    return (a - b) / abs(b) if b and abs(b) > 1e-15 else np.nan


def scores(tag):
    f = f'{V2}/results/scores/score_series_{tag}.csv'
    return pd.read_csv(f).set_index('site') if os.path.isfile(f) else None


def compare(case, ref, tag):
    done = lambda t, s: os.path.isfile(f'{R}/cases/{case}/sites_{t}/{s}/.complete')
    sites = sorted(os.path.basename(os.path.dirname(p)) for p in
                   glob.glob(f'{R}/cases/{case}/sites_{ref}/*/.complete'))
    sites = [s for s in sites if done(tag, s)]
    sr, st = scores(ref), scores(tag)
    rows = []
    for s in sites:
        e0, e1 = annual(case, ref, s, CH4, 'tracer_'), annual(case, tag, s, CH4, 'tracer_')
        yrs = sorted(set(e0) & set(e1))
        if not yrs:
            continue
        d = {'site': s, 'dE_y1': rel(e1[yrs[0]], e0[yrs[0]]), 'dE_last': rel(e1[yrs[-1]], e0[yrs[-1]]),
             'dE_all': rel(np.mean([e1[y] for y in yrs]), np.mean([e0[y] for y in yrs]))}
        for key, var in (('dSOM', 'f_totsomc'), ('dLIT', 'f_totlitc')):
            p0, p1 = annual(case, ref, s, var, ''), annual(case, tag, s, var, '')
            y = sorted(set(p0) & set(p1))
            d[key] = rel(p1[y[0]], p0[y[0]]) if y else np.nan
        if sr is not None and st is not None and s in sr.index and s in st.index:
            for k in ('KGEln', 'beta', 'alpha', 'r'):
                d['d' + k] = st.at[s, k] - sr.at[s, k]
        rows.append(d)
    return pd.DataFrame(rows).set_index('site')


def main():
    case, ref, tags = sys.argv[1], sys.argv[2], sys.argv[3:]
    out = []
    for tag in tags:
        df = compare(case, ref, tag)
        pct = [c for c in df.columns if c.startswith(('dE', 'dSOM', 'dLIT'))]
        show = df.copy()
        show[pct] = show[pct] * 100
        out.append(f'== {tag} against {ref} ({case}); dE, dSOM, dLIT in %; {len(df)} towers')
        out.append(show.round(2).to_string())
        med = show.abs().median()
        worst = show.abs().max()
        out.append('median |d|: ' + '  '.join(f'{c} {med[c]:.2f}' for c in show.columns))
        out.append('max    |d|: ' + '  '.join(f'{c} {worst[c]:.2f}' for c in show.columns))
        ns = df[~df.index.isin(SALT)]
        if 'dKGEln' in ns:
            out.append(f'20 non-saline, median score change: KGE_ln {ns.dKGEln.median():+.3f}  '
                       f'beta {ns.dbeta.median():+.3f}  alpha {ns.dalpha.median():+.3f}  r {ns.dr.median():+.3f}')
        out.append('')
    text = '\n'.join(out)
    print(text)
    os.makedirs(f'{V2}/results', exist_ok=True)
    with open(f'{V2}/results/eq_{ref}__{"_".join(tags)}.txt', 'w') as f:
        f.write(text + '\n')


if __name__ == '__main__':
    main()
