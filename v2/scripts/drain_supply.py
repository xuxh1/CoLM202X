#!/usr/bin/env python
"""Annual lateral drainage (positive f_rsub) and floor supply (negative f_rsub)
of the dynamic wetland per tower for two site tags, from monthly means
(mm/yr; a month holding both shows only its net). Usage: drain_supply.py
<case_a> <tag_a> <case_b> <tag_b>"""
import csv
import glob
import sys

import netCDF4 as nc
import numpy as np

R = '/share/home/dq076/mode/Methane'
SITES = [r['ID'] for r in csv.DictReader(open(f'{R}/scripts/sites/LIST_sites_wtd23.csv'))]


def ann(case, tag, s):
    pos, neg, n = 0.0, 0.0, 0
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as h:
            r = np.ma.filled(h['f_rsub'][:], np.nan).ravel()
            p = np.ma.filled(h['f_xy_rain'][:], np.nan).ravel() + np.ma.filled(h['f_xy_snow'][:], np.nan).ravel()
        if r.size != 12:
            continue
        sec = 86400 * 30.44
        pos += np.nansum(np.clip(r, 0, None)) * sec
        neg += np.nansum(np.clip(r, None, 0)) * sec
        n += 1
    return pos / n, -neg / n


def main():
    ca, ta, cb, tb = sys.argv[1:5]
    print(f'# mm/yr; drain = positive f_rsub, supply = negative f_rsub; {ta} vs {tb}')
    print(f'{"site":8s} {"drainA":>8s} {"supplyA":>8s} {"drainB":>8s} {"supplyB":>8s}')
    for s in SITES:
        a, b = ann(ca, ta, s), ann(cb, tb, s)
        print(f'{s:8s} {a[0]:8.0f} {a[1]:8.0f} {b[0]:8.0f} {b[1]:8.0f}')


if __name__ == '__main__':
    main()
