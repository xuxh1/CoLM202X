"""Score the rice towers of one site tag against the paddy-only CH4 flux (V-25).

A rice tower's footprint is paddy, but a single-point run on CROPLAND takes
every crop type of the coarse CFT map at the tower (US-HRA: 77% rice, IT-Cas
30%), so the column-mean flux mixes paddy and upland crops. This scores the
paddy-area intensive flux f_methane_surf_flux_rice_intensive against the
tower with score_series' monthly pairing and KGE terms, and adds the paddy
budget: emission over production (E/P) and the ebullition / aerenchyma /
diffusion shares of the channel sum, all over the tower's observed years.
Usage: rice_sites.py <case dir> <tag> [site ...]"""
import glob
import os
import re
import sys

import numpy as np
import netCDF4 as nc

sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/sites')
import score_series as S                                              # noqa: E402

SITES = ['IT-Cas', 'JP-Mse', 'KR-CRK', 'PH-RiF', 'US-HRA', 'US-HRC', 'US-Twt']
BUDGET = ('f_methane_prod_tot_rice', 'f_methane_oxid_tot_rice', 'f_methane_surf_ebul_rice',
          'f_methane_surf_aere_rice', 'f_methane_surf_diff_rice', 'f_methane_area_rice')


def monthly(case, tag, sid, var):
    """(year, month) -> monthly value of one tracer variable."""
    out = {}
    for p in sorted(glob.glob(f'{case}/sites_{tag}/{sid}/history/*_hist_*tracer*.nc')):
        m = re.search(r'_(\d{4})[_.]', os.path.basename(p))
        if not m:
            continue
        with nc.Dataset(p) as d:
            if var not in d.variables:
                continue
            v = np.ma.filled(d.variables[var][:].astype(float), np.nan)
            # single-point history is (time, patch): take the rice patch, i.e. the
            # column whose f_methane_area_rice is ever positive (never flatten)
            v = v.reshape(v.shape[0], -1)
            col = 0
            if v.shape[1] > 1 and 'f_methane_area_rice' in d.variables:
                a = np.ma.filled(d.variables['f_methane_area_rice'][:].astype(float), 0).reshape(v.shape[0], -1)
                col = int(np.argmax(a.max(0)))
            v = v[:, col]
        if v.size == 12:
            for i in range(12):
                out[(int(m.group(1)), i + 1)] = float(v[i])
    return out


def main():
    case, tag = sys.argv[1:3]
    sites = sys.argv[3:] or SITES
    print(f'# {tag} ({case}); rice towers vs paddy-intensive CH4 flux, monthly pairs')
    print('site     n  KGEln   beta  alpha      r | obs_mean mod_mean mg/m2/d | '
          'E/P  ebul aere diff | paddy_frac')
    for sid in sites:
        if not os.path.isdir(f'{case}/sites_{tag}/{sid}'):
            print(f'{sid:7s} no run')
            continue
        obs = S.obs_monthly(sid)
        mod = {k: v * S.MOD_TO_MGM2D for k, v in
               monthly(case, tag, sid, 'f_methane_surf_flux_rice_intensive').items()}
        keys = sorted(k for k in obs if k in mod and np.isfinite(obs[k]) and np.isfinite(mod[k]))
        if len(keys) < 3:
            print(f'{sid:7s} {len(keys):2d} pairs')
            continue
        o = np.array([obs[k] for k in keys])
        m = np.array([mod[k] for k in keys])
        _, kln, beta, alpha, r = S.kge_terms(o, m)
        years = {k[0] for k in keys}
        b = {v: monthly(case, tag, sid, v) for v in BUDGET}
        tot = {v: np.nansum([x for k, x in b[v].items() if k[0] in years]) for v in BUDGET[:5]}
        chan = tot['f_methane_surf_ebul_rice'] + tot['f_methane_surf_aere_rice'] + tot['f_methane_surf_diff_rice']
        ep = chan / tot['f_methane_prod_tot_rice'] if tot['f_methane_prod_tot_rice'] > 0 else np.nan
        sh = [tot[v] / chan if chan > 0 else np.nan for v in BUDGET[2:5]]
        frac = np.nanmax(list(b['f_methane_area_rice'].values()) or [np.nan])
        print(f'{sid:7s} {len(keys):2d} {kln:6.2f} {beta:6.2f} {alpha:6.2f} {r:6.2f} | '
              f'{o.mean():8.1f} {m.mean():8.1f}        | {ep:4.2f} {sh[0]:5.2f} {sh[1]:4.2f} {sh[2]:4.2f} | {frac:5.2f}')


if __name__ == '__main__':
    main()
