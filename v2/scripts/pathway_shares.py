#!/usr/bin/env python
"""Internal CH4 partition of the model: fraction of production oxidized and the
shares of emission via aerenchyma, ebullition and diffusion (annual sums).
sites <case> <tag>        : per tower (inundated-fraction weighted sat/unsat
                            subcolumns), with class medians over the 20
                            non-saline towers of LIST_sites_wtd23.
global <ver> <name> <year>: wetland tile and floodplain (soil tile) by band.
Usage: pathway_shares.py sites <case under cases/> <tag> | global <ver> <name> <year>"""
import csv
import glob
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
CLS = {1: 'permafrost', 3: 'bog', 4: 'fen', 5: 'marsh', 6: 'salt_marsh', 7: 'trop_swamp'}
Q = ('prod', 'oxid', 'aere', 'ebul', 'diff')


def site_sums(case, tag, s):
    acc = dict.fromkeys(Q, 0.0)
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_tracer_2*.nc')):
        with nc.Dataset(f) as e:
            F = np.ma.filled(e['f_methane_finundated'][:], np.nan).ravel()
            if F.size != 12:
                continue
            for k, a, b in (('prod', 'f_methane_prod_tot_sat', 'f_methane_prod_tot_unsat'),
                            ('oxid', 'f_methane_oxid_tot_sat', 'f_methane_oxid_tot_unsat'),
                            ('aere', 'f_methane_surf_aere_sat', 'f_methane_surf_aere_unsat'),
                            ('ebul', 'f_methane_surf_ebul_sat', 'f_methane_surf_ebul_unsat'),
                            ('diff', 'f_methane_surf_diff_sat', 'f_methane_surf_diff_unsat')):
                A = np.ma.filled(e[a][:], np.nan).ravel(); Bv = np.ma.filled(e[b][:], np.nan).ravel()
                acc[k] += float(np.nansum(F * A + (1 - F) * Bv))
    return acc


def shares(a):
    em = a['aere'] + a['ebul'] + a['diff']
    return {'oxid/prod': a['oxid'] / a['prod'] if a['prod'] > 0 else np.nan,
            'aere/emis': a['aere'] / em if em > 0 else np.nan,
            'ebul/emis': a['ebul'] / em if em > 0 else np.nan,
            'diff/emis': a['diff'] / em if em > 0 else np.nan}


def sites(case, tag):
    rows = []
    for r in csv.DictReader(open(f'{R}/scripts/sites/LIST_sites_wtd23.csv')):
        d = shares(site_sums(case, tag, r['ID']))
        d.update(site=r['ID'], cls=CLS[int(r['SITE_wetland_class'])])
        rows.append(d)
    t = pd.DataFrame(rows).set_index('site')
    print(f'# {tag}: annual shares per tower')
    print(t.round(2).to_string())
    print('# class medians (salt marsh excluded from the 20-tower median)')
    print(t.groupby('cls').median(numeric_only=True).round(2).to_string())
    print('20 non-saline median:', t[t.cls != 'salt_marsh'].median(numeric_only=True).round(2).to_dict())


def global_(ver, name, y):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import v2budget as VB
    import xarray as xr
    B = VB.B
    c = B.Case(ver, name)
    A_land, _, _ = VB.load_areas(c, y)
    with xr.open_dataset(f'{c.dir}/history/{name}_hist_tracer_{y}.nc') as tr:
        w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')
        tot = lambda v: (tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4
        for tile in ('wetland', 'soil'):
            d = {k: tot(f'f_methane_{"prod_tot" if k == "prod" else "oxid_tot" if k == "oxid" else "surf_" + k}_{tile}')
                 for k in Q}
            for b, lo, hi in (('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90), ('all', -90, 90)):
                m = (d['prod'].lat >= lo) & (d['prod'].lat < hi)
                a = {k: float(v.where(m).sum()) for k, v in d.items()}
                sh = shares(a)
                print(f'{tile:8s} {b:8s} prod {a["prod"]:7.1f} oxid {a["oxid"]:6.1f} emis {a["aere"] + a["ebul"] + a["diff"]:6.1f} Tg | '
                      + '  '.join(f'{k} {v:4.2f}' for k, v in sh.items()))


if __name__ == '__main__':
    if sys.argv[1] == 'sites':
        sites(sys.argv[2], sys.argv[3])
    else:
        global_(sys.argv[2], sys.argv[3], int(sys.argv[4]))
