#!/usr/bin/env python
"""Bitwise comparison of two site tags: every variable of every history file
(main and tracer) of each site common to both. Usage:
    compare_hist.py <case_a>/sites_<tag_a> <case_b>/sites_<tag_b>
Prints per site the number of variables compared and those that differ."""
import glob
import os
import sys

import netCDF4 as nc
import numpy as np


def main():
    a, b = sys.argv[1], sys.argv[2]
    sites = sorted(set(os.listdir(a)) & set(os.listdir(b)) - {'_conf'})
    bad = 0
    for s in sites:
        fa = sorted(glob.glob(f'{a}/{s}/history/*.nc'))
        n, diff = 0, []
        for f in fa:
            g = f.replace(a, b, 1)
            if not os.path.exists(g):
                diff.append(f'missing {os.path.basename(g)}')
                continue
            with nc.Dataset(f) as x, nc.Dataset(g) as y:
                for v in x.variables:
                    if v not in y.variables:
                        diff.append(f'{os.path.basename(f)}:{v} missing')
                        continue
                    p, q = np.ma.filled(x[v][:], np.nan), np.ma.filled(y[v][:], np.nan)
                    n += 1
                    if p.shape != q.shape or not np.array_equal(p, q, equal_nan=True):
                        diff.append(f'{os.path.basename(f)}:{v}')
        bad += len(diff)
        print(f'{s}: {n} variables compared, {len(diff)} differ' + (': ' + ', '.join(diff[:8]) if diff else ''))
    print('IDENTICAL' if bad == 0 else f'DIFFERENT ({bad})')


if __name__ == '__main__':
    main()
