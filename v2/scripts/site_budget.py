"""Carbon and CH4 budget of each tower in a site tag against the tower (main-set
diagnosis): over the months both have, the model's heterotrophic respiration
HR, CH4 production P, oxidation O and emission E, the saturated share of the
potential decomposition, and E/HR; the tower's CH4 emission and ecosystem
respiration Reco (FCH4_f, Resp_DT, gap-filled, months with at least 80% of the
steps) and E/(Reco/2), taking half of Reco as heterotrophic (assumption).
Fluxes in g C m-2 yr-1.
Usage: site_budget.py <case under cases/> <tag>"""
import glob
import re
import sys
from collections import defaultdict
import numpy as np
import netCDF4 as nc

C = '/share/home/dq076/mode/Methane/cases/'
OBS = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Observation/'
SPY = 365.25 * 86400
MV = {'HR': ('', 'f_hr', 1.), 'P': ('_tracer', 'f_methane_prod_tot', 12.),
      'O': ('_tracer', 'f_methane_oxid_tot', 12.), 'E': ('_tracer', 'f_methane_surf_flux_tot_active', 12.),
      'Dsat': ('_tracer', 'f_co2_decomp_tot_sat', 12.), 'Dunsat': ('_tracer', 'f_co2_decomp_tot_unsat', 12.)}


def obs_monthly(s):
    acc = {'E': defaultdict(list), 'R': defaultdict(list)}
    for p in glob.glob(f'{OBS}{s}_*_Flux.nc'):
        with nc.Dataset(p) as d:
            t = nc.num2date(d['time'][:], d['time'].units)
            keys = [(x.year, x.month) for x in t]
            for k, v, fac in (('E', 'FCH4_f', 1e-9 * 12), ('R', 'Resp_DT', 1e-6 * 12)):
                if v not in d.variables:
                    continue
                a = np.ma.filled(d[v][:].astype(float), np.nan).ravel() * fac * SPY
                for kk, vv in zip(keys, a):
                    acc[k][kk].append(vv)
    return {k: {kk: np.nanmean(v) for kk, v in m.items() if np.isfinite(v).mean() >= 0.8} for k, m in acc.items()}


def mod_monthly(case, tag, s):
    out = {k: {} for k in MV}
    for f in sorted(glob.glob(f'{C}{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        y = int(re.search(r'_(\d{4})\.nc$', f).group(1))
        for k, (suf, v, fac) in MV.items():
            try:
                with nc.Dataset(f.replace('_hist_', f'_hist{suf}_')) as d:
                    a = np.ma.filled(d[v][:].astype(float), np.nan).reshape(d[v].shape[0], -1)[:, 0] * fac * SPY
                    # a tower's first year may start mid-year (ID-Pag: June);
                    # records carry mid-month time stamps
                    mon = [x.month for x in nc.num2date(d['time'][:], d['time'].units)]
            except (OSError, IndexError):
                continue
            for i in range(min(len(mon), a.size)):
                out[k][(y, mon[i])] = a[i]
    return out


def main():
    case, tag = sys.argv[1:3]
    sites = sorted(p.split('/')[-1] for p in glob.glob(f'{C}{case}/sites_{tag}/*') if '_conf' not in p)
    print(f'# {tag} ({case}); g C m-2 yr-1 over months with tower CH4 (Reco columns: months with both)')
    print('site     n |   HR     P     O     E  sat%  P/HR%  E/HR% | E_obs  Reco  E_obs/(Reco/2)% | E/E_obs')
    for s in sites:
        o, m = obs_monthly(s), mod_monthly(case, tag, s)
        ks = [k for k in o['E'] if all(k in m[x] for x in MV)]
        if len(ks) < 6:
            print(f'{s:7s} {len(ks):2d}')
            continue
        a = {x: np.mean([m[x][k] for k in ks]) for x in MV}
        eo = np.mean([o['E'][k] for k in ks])
        kr = [k for k in ks if k in o['R']]
        ro = np.mean([o['R'][k] for k in kr]) if len(kr) >= 6 else np.nan
        er = np.mean([m['E'][k] for k in kr]) / np.mean([m['HR'][k] for k in kr]) if len(kr) >= 6 else np.nan
        sat = a['Dsat'] / (a['Dsat'] + a['Dunsat']) if a['Dsat'] + a['Dunsat'] > 0 else np.nan
        eor = np.mean([o['E'][k] for k in kr]) / (ro / 2) if len(kr) >= 6 else np.nan
        print(f"{s:7s} {len(ks):2d} | {a['HR']:5.0f} {a['P']:5.1f} {a['O']:5.1f} {a['E']:5.1f} {100 * sat:4.0f} "
              f"{100 * a['P'] / a['HR']:6.1f} {100 * a['E'] / a['HR']:6.1f} | {eo:5.1f} {ro:5.0f} {100 * eor:6.1f} | "
              f"{a['E'] / eo:5.2f}  (model E/HR on Reco months {100 * er:.1f}%)")


if __name__ == '__main__':
    main()
