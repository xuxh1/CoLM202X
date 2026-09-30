#!/usr/bin/env python
"""Winter storage and thaw release at cold towers: monthly climatology of the
column CH4 stock (saturated and unsaturated subcolumns weighted by the
inundated fraction), production, oxidation, and surface emission by pathway
(aerenchyma, ebullition, diffusion), with top-layer soil temperature, snow
depth, inundated fraction and observed CH4. Units mg CH4 m-2 d-1, stock
mg CH4 m-2. Usage: winter_budget.py <case under cases/> <tag> <site>..."""
import glob
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
sys.path.insert(0, f'{R}/scripts/sites')
import score_series as SS                                           # noqa: E402
MG = 16.04e3                                                        # mol -> mg CH4


def clim(case, tag, s):
    rows = []
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as h, nc.Dataset(f.replace('_hist_', '_hist_tracer_')) as e:
            g = lambda d, v: np.ma.filled(d[v][:], np.nan).astype(float) if v in d.variables else None
            flux = g(e, 'f_methane_surf_flux_tot_active')
            if flux is None or flux.size != 12:
                continue
            F = g(e, 'f_methane_finundated').ravel()
            wsum = lambda a, b: F * a.ravel() + (1 - F) * b.ravel()
            T = g(h, 'f_t_soisno').reshape(12, -1)[:, -10:][:, 0] - 273.15
            snow = g(h, 'f_snowdp')
            for i in range(12):
                r = {'m': i + 1, 'Ttop': T[i], 'F': F[i], 'emis': flux.ravel()[i] * SS.MOD_TO_MGM2D,
                     'snow_m': snow.ravel()[i] if snow is not None else np.nan}
                for k, a, b in (('stock', 'f_totcol_methane_sat', 'f_totcol_methane_unsat'),
                                ('prod', 'f_methane_prod_tot_sat', 'f_methane_prod_tot_unsat'),
                                ('oxid', 'f_methane_oxid_tot_sat', 'f_methane_oxid_tot_unsat'),
                                ('aere', 'f_methane_surf_aere_sat', 'f_methane_surf_aere_unsat'),
                                ('ebul', 'f_methane_surf_ebul_sat', 'f_methane_surf_ebul_unsat'),
                                ('diff', 'f_methane_surf_diff_sat', 'f_methane_surf_diff_unsat')):
                    A, Bv = g(e, a), g(e, b)
                    if A is None or Bv is None:
                        r[k] = np.nan
                        continue
                    v = wsum(A, Bv)[i]
                    r[k] = v * MG if k == 'stock' else v * SS.MOD_TO_MGM2D
                rows.append(r)
    return pd.DataFrame(rows).groupby('m').mean()


def main():
    case, tag = sys.argv[1], sys.argv[2]
    for s in sys.argv[3:]:
        c = clim(case, tag, s)
        c['obs'] = pd.Series(SS.obs_monthly(s)).groupby(level=1).mean()
        print(f'== {s} ({tag})')
        print(c[['obs', 'emis', 'prod', 'oxid', 'aere', 'ebul', 'diff', 'stock', 'F', 'Ttop', 'snow_m']].round(2).T.to_string())


if __name__ == '__main__':
    main()
