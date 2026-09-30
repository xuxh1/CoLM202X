"""Cold-season CH4 at the towers (Q-27): monthly climatology of the tower flux
and of the model's column production P, oxidation O, surface flux by channel
(diffusion, ebullition, aerenchyma), column CH4 store, snow depth and the
temperature and ice share of soil layers 1, 3 and 5 (nodes about 0.7, 6 and
21 cm), to tell whether autumn and winter emission stop because production
stops or because the CH4 cannot leave the column.
Usage: cold_season.py <case under cases/> <tag> <site> [site ...]"""
import glob
import sys
from collections import defaultdict
import numpy as np
import netCDF4 as nc

sys.path.insert(0, '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/scripts')
import score_series_v2 as S                                          # noqa: E402

C = '/share/home/dq076/mode/Methane/cases/'
K = 16.043e3 * 86400                     # mol m-2 s-1 -> mg CH4 m-2 d-1
TV = {'P': 'f_methane_prod_tot', 'O': 'f_methane_oxid_tot', 'E': 'f_methane_surf_flux_tot_active',
      'diff': 'f_methane_surf_diff_wetland', 'ebul': 'f_methane_surf_ebul_wetland',
      'aere': 'f_methane_surf_aere_wetland', 'store': 'f_totcol_methane'}


def main():
    case, tag = sys.argv[1:3]
    for sid in sys.argv[3:]:
        acc = defaultdict(lambda: defaultdict(list))
        for f in sorted(glob.glob(f'{C}{case}/sites_{tag}/{sid}/history/{sid}_hist_2*.nc')):
            with nc.Dataset(f) as d, nc.Dataset(f.replace('_hist_', '_hist_tracer_')) as t:
                n = d['f_zwt'].shape[0]
                if n < 12:
                    continue
                g = lambda ds, v: np.ma.filled(ds[v][:].astype(float), np.nan).reshape(n, -1)[:, 0]
                for k, v in TV.items():
                    if v in t.variables:
                        a = g(t, v) * (K if k != 'store' else 16.043)          # store in g CH4 m-2
                        for m in range(12):
                            acc[k][m].append(a[m])
                ts = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)[:, 0, :]
                wl = np.ma.filled(d['f_wliq_soisno'][:].astype(float), np.nan)[:, 0, :]
                wi = np.ma.filled(d['f_wice_soisno'][:].astype(float), np.nan)[:, 0, :]
                sd = g(d, 'f_snowdp')
                for m in range(12):
                    acc['snow'][m].append(sd[m])
                    for L, j in (('1', 5), ('3', 7), ('5', 9)):
                        acc['T' + L][m].append(ts[m, j] - 273.15)
                        acc['ice' + L][m].append(wi[m, j] / (wi[m, j] + wl[m, j]) if wi[m, j] + wl[m, j] > 0 else np.nan)
        obs = S.obs_monthly(sid)
        om = defaultdict(list)
        for (y, m), v in obs.items():
            if np.isfinite(v):
                om[m - 1].append(v)
        print(f'== {sid} ({tag})  month ' + ' '.join(f'{m:6d}' for m in range(1, 13)))
        print('  obs E            ' + ' '.join(f'{np.mean(om[m]):6.1f}' if om[m] else '     .' for m in range(12)))
        for k, fmt in (('E', '6.1f'), ('P', '6.1f'), ('O', '6.1f'), ('diff', '6.1f'), ('ebul', '6.1f'), ('aere', '6.1f'),
                       ('store', '6.2f'), ('snow', '6.2f'), ('T1', '6.1f'), ('T3', '6.1f'), ('T5', '6.1f'),
                       ('ice1', '6.2f'), ('ice3', '6.2f'), ('ice5', '6.2f')):
            if acc[k]:
                print(f'  model {k:10s} ' + ' '.join(f'{np.nanmean(acc[k][m]):{fmt}}' if acc[k][m] else '     .' for m in range(12)))


if __name__ == '__main__':
    main()
