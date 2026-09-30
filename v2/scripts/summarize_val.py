#!/usr/bin/env python
"""Score one V2 validation tag (the 23 wetland towers outside the main set,
v2/config/LIST_sites_val23.csv) and write the summary into v2/results/.
Same steps as summarize_tag.py: score_series.py on the validation list, CH4
medians of KGE_ln, beta, alpha, r over all towers and per tower-label class;
water level (ledger D-2) where observed depths exist.
Usage: summarize_val.py <case under cases/> <tag>"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize_tag import OB, PY, R, V2, water_level                # noqa: E402

LIST = f'{V2}/config/LIST_sites_val23.csv'


def main():
    case, tag = sys.argv[1], sys.argv[2]
    subprocess.run([PY, '-W', 'ignore', f'{V2}/scripts/score_series_v2.py', tag, f'{R}/cases/{case}',
                    '--sitelist', LIST], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.makedirs(f'{V2}/results/scores', exist_ok=True)
    for ext in ('txt', 'csv'):
        shutil.move(f'{R}/outputs/analysis/score_series_{tag}.{ext}', f'{V2}/results/scores/score_series_{tag}.{ext}')
    s = pd.read_csv(f'{V2}/results/scores/score_series_{tag}.csv')
    lab = pd.read_csv(LIST).set_index('ID').SITE_CLASSIFICATION
    s['label'] = s.site.map(lab)
    four = lambda d: ' / '.join(f'{d[k].median():.2f}' for k in ('KGEln', 'beta', 'alpha', 'r'))
    L = [f'# {tag} ({case}); validation towers; CH4 medians KGE_ln / beta / alpha / r',
         f'all {len(s)}: {four(s)}']
    for c, d in s.groupby('label'):
        L.append(f'  {c:11s} n={len(d):2d}  {four(d)}')
    wl = []
    for site in s.site:
        m = water_level(case, tag, site)
        if m.empty:
            continue
        j = m.merge(OB[OB.site == site][['year', 'month', 'wtd_obs']], on=['year', 'month']).dropna()
        if len(j) < 6:
            continue
        r = j.lev.corr(j.wtd_obs) if j.lev.std() > 1e-6 else np.nan
        wl.append({'site': site, 'obs': j.wtd_obs.median(), 'mod': j.lev.median(),
                   'err': abs(j.lev.median() - j.wtd_obs.median()), 'r': r})
    w = pd.DataFrame(wl)
    if len(w):
        L.append(f'water level (D-2), {len(w)} towers: median |median error| {w.err.median():.2f} m, '
                 f'median r {w.r.median():.2f}')
        L.append(w.round(2).to_string(index=False))
    else:
        L.append('water level (D-2): no tower with observed depths')
    L.append('')
    L.append(s[['site', 'label', 'KGEln', 'beta', 'alpha', 'r']].round(2).to_string(index=False))
    open(f'{V2}/results/{tag}_summary.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:10]))


if __name__ == '__main__':
    main()
