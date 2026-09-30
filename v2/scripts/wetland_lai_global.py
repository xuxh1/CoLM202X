#!/usr/bin/env python
"""Peak monthly LAI (2005 surface data) of wetland patches (IGBP class 11) by
20-degree latitude block, GLWD surface data (V-3/V-6) against the paper run's
surface data, with the soil patch (class 1, all natural vegetation) of the same blocks for reference.
Patch order in LAI_patchesMM_<block>.nc follows landpatch_<block>.nc."""
import glob
import os

import netCDF4 as nc
import numpy as np

CASES = {'GLWD (V-3/V-6)': '/share/home/dq076/mode/Methane/cases/paper_v2/v260924c/g2_v2/landdata',
         'paper run': '/share/home/dq076/mode/Methane/cases/v260910c/g2_paper_fp2/landdata'}
BANDS = {'n70': '70-90N', 'n50': '50-70N', 'n30': '30-50N', 'n10': '10-30N', 's10': '10S-10N', 's30': '30-10S', 's50': '50-30S'}
GROUP = {'wetland (11)': [11], 'soil patch (1)': [1]}
for lab, L in CASES.items():
    L = os.path.realpath(L)
    out = {}
    for f in sorted(glob.glob(f'{L}/landpatch/2005/landpatch_*.nc')):
        blk = f.split('landpatch_')[-1][:-3]
        band = blk.split('_')[1]
        if band not in BANDS:
            continue
        typ = np.asarray(nc.Dataset(f)['settyp'][:])
        lai = []
        for m in range(1, 13):
            g = f'{L}/LAI/2005/LAI_patches{m:02d}_{blk}.nc'
            if not os.path.exists(g):
                lai = None
                break
            with nc.Dataset(g) as h:
                lai.append(np.ma.filled(h['LAI_patches'][:], np.nan).astype(float))
        if lai is None or len(lai[0]) != len(typ):
            continue
        peak = np.nanmax(np.array(lai), axis=0)
        for gname, cl in GROUP.items():
            out.setdefault((band, gname), []).extend(peak[np.isin(typ, cl)].tolist())
    print(f'== {lab}: peak monthly LAI of patches, median [p25-p75] (n)')
    for band, bl in BANDS.items():
        row = []
        for gname in GROUP:
            v = np.array(out.get((band, gname), []))
            v = v[np.isfinite(v) & (v < 50)]
            row.append(f'{gname} {np.median(v):4.2f} [{np.percentile(v, 25):4.2f}-{np.percentile(v, 75):4.2f}] ({len(v)})'
                       if len(v) else f'{gname} none')
        print(f'  {bl:8s} ' + ' | '.join(row))
