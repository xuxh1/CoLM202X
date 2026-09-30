#!/usr/bin/env python
"""Score a V2 site tag on all 77 towers of scripts/sites/LIST_sites_77.csv
(model-optimisation loop, user 2026-09-27: every tower, not only the 23 with
a water table).

Observation gate and months: score_series_v2 (ledger D-17). Model monthly flux:
  lake towers      f_methane_surf_flux_tot (land-area mean incl. the lake; the
                   tower cell is all lake), the active-area mean excludes lakes
  other towers     f_methane_surf_flux_tot_active
  several patches  weighted by pctcrop of the tower's landdata/srfdata.nc (crop
                   towers without a single-crop patch); equal weights otherwise
Per tower: n months, KGE_ln, beta, alpha, r, observed and model mean
(mg CH4 m-2 d-1). Classes are the FLUXNET-CH4 labels of the list, with US-Uaf
reported as 'permafrost' and the brackish DE-Hte with the salt marshes (as
all_sites.py and summarize_tag.py). Medians of the four numbers per
class and over
  wetland  = bog, fen, marsh, swamp, wet tundra, drained, permafrost (the
             former 38-tower set plus the new drained towers)
  all      = wetland + rice + lake
Upland towers (net fluxes near zero or negative) are listed with the observed
and model means and r, and summarised by the median of (model - obs) and the
share with the right sign; salt marshes are listed apart (no sulphate
suppression in the model).
Writes v2/results/scores77_<tag>.txt and .csv.
Usage: score77.py <case under cases/> <tag>"""
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_series_v2 as V                                         # noqa: E402

S = V.S
R = '/share/home/dq076/mode/Methane'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
LIST = f'{R}/scripts/sites/LIST_sites_77.csv'
WET = ('bog', 'fen', 'marsh', 'swamp', 'wet tundra', 'drained', 'permafrost')


def classes():
    c = {}
    for row in csv.DictReader(open(LIST)):
        c[row['ID']] = row['SITE_CLASSIFICATION'].strip().lower()
    c['US-Uaf'] = 'permafrost'
    c['DE-Hte'] = 'salt marsh'          # brackish coastal fen, V2 keeps it with the salt marshes
    return c


def weights(case, tag, sid, npatch):
    if npatch == 1:
        return np.ones(1)
    f = f'{R}/cases/{case}/sites_{tag}/{sid}/landdata/srfdata.nc'
    if os.path.exists(f):
        with nc.Dataset(f) as d:
            if 'pctcrop' in d.variables:
                w = np.ma.filled(d['pctcrop'][:], 0.0).astype(float).ravel()
                if w.size == npatch and w.sum() > 0:
                    return w / w.sum()
    return np.ones(npatch) / npatch


def mod_monthly(case, tag, sid, lake):
    var = 'f_methane_surf_flux_tot' if lake else S.CH4
    out = {}
    for p in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{sid}/history/*_hist_*tracer*.nc')):
        m = re.search(r'_(\d{4})[_.]', os.path.basename(p))
        if not m:
            continue
        with nc.Dataset(p) as d:
            if var not in d.variables:
                continue
            v = np.ma.filled(d.variables[var][:], np.nan).astype(float)
        if v.ndim == 1:
            v = v[:, None]
        if v.shape[0] != 12:
            continue
        v = np.nansum(v * weights(case, tag, sid, v.shape[1])[None, :], 1) * S.MOD_TO_MGM2D
        for i in range(12):
            out[(int(m.group(1)), i + 1)] = float(v[i])
    return out


def footprint():
    """Wetland share of the flux footprint of towers whose footprint is mostly
    not wetland (v2/config/FOOTPRINT_WETLAND.tsv, user 2026-09-29, Q-50); the
    modelled wetland flux is scaled by it before scoring."""
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'FOOTPRINT_WETLAND.tsv')
    out = {}
    if os.path.exists(f):
        for line in open(f):
            if line.strip() and not line.startswith('#'):
                p = line.split('\t')
                out[p[0].strip()] = float(p[1])
    return out


def med(rows, k):
    a = np.array([r[k] for r in rows], float)
    a = a[np.isfinite(a)]
    return np.median(a) if a.size else np.nan


def main():
    case, tag = sys.argv[1:3]
    cls = classes()
    fp = footprint()
    rows, skip = [], []
    for sid in sorted(cls):
        if S.obs_clim_months(sid) < 8:
            skip.append(f'{sid} (<8 obs months)')
            continue
        om = V.obs_monthly(sid)
        mm = mod_monthly(case, tag, sid, cls[sid] == 'lake')
        if sid in fp:
            mm = {k: v * fp[sid] for k, v in mm.items()}
        keys = sorted(set(om) & set(mm))
        if len(keys) < 8:
            skip.append(f'{sid} ({len(keys)} common months)')
            continue
        o = np.array([om[k] for k in keys])
        m = np.array([mm[k] for k in keys])
        _, kln, beta, alpha, r = S.kge_terms(o, m)
        rows.append(dict(site=sid, cls=cls[sid], n=len(keys), kln=kln, beta=beta, alpha=alpha,
                         r=r, obs=o.mean(), mod=m.mean()))
    L = [f'# {case} {tag}: all 77 towers; KGE_ln / beta / alpha / r; means mg CH4 m-2 d-1',
         f'{"site":8s} {"class":11s} {"n":>3s} {"KGEln":>6s} {"beta":>6s} {"alpha":>6s} {"r":>5s} {"obs":>7s} {"model":>7s}']
    order = WET + ('rice', 'lake', 'salt marsh', 'upland')
    for c in order:
        for x in [x for x in rows if x['cls'] == c]:
            L.append(f'{x["site"]:8s} {x["cls"]:11s} {x["n"]:3d} {x["kln"]:6.2f} {x["beta"]:6.2f} '
                     f'{x["alpha"]:6.2f} {x["r"]:5.2f} {x["obs"]:7.2f} {x["mod"]:7.2f}')
    L.append('# unscored: ' + ', '.join(skip))
    L.append('')
    L.append(f'{"group":12s} {"n":>3s} {"KGEln":>6s} {"beta":>6s} {"alpha":>6s} {"r":>5s}')

    def line(name, sel):
        if sel:
            L.append(f'{name:12s} {len(sel):3d} ' + ' '.join(f'{med(sel, k):6.2f}' for k in ('kln', 'beta', 'alpha', 'r')))

    wet = [x for x in rows if x['cls'] in WET]
    line('all', wet + [x for x in rows if x['cls'] in ('rice', 'lake')])
    line('wetland', wet)
    for c in order[:-1]:
        line(c, [x for x in rows if x['cls'] == c])
    up = [x for x in rows if x['cls'] == 'upland']
    if up:
        d = np.array([x['mod'] - x['obs'] for x in up])
        sgn = np.mean([np.sign(x['mod']) == np.sign(x['obs']) for x in up])
        L.append(f'upland       {len(up):3d} median(model-obs) {np.median(d):6.2f} mg m-2 d-1; '
                 f'right sign {sgn:.2f}; median r {med(up, "r"):5.2f}')
    txt = '\n'.join(L)
    print(txt)
    open(f'{OUT}/scores77_{tag}.txt', 'w').write(txt + '\n')
    with open(f'{OUT}/scores77_{tag}.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
