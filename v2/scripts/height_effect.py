"""What the tower measurement heights changed at the sites (V-29): per tower,
the annual means of friction velocity, sensible and latent heat, ground heat,
soil temperature at about 10 cm and 30 cm, snow depth, water table, and the
CH4 flux, 20 m run (v2_k4p) against the tower-height run (v2_k4ph)."""
import glob, sys
import numpy as np, netCDF4 as nc
B = '/share/home/dq076/mode/Methane/cases/paper_v2/v260925u/sp_v2_c17c/'
V = ['f_ustar', 'f_fsena', 'f_lfevpa', 'f_fgrnd', 'f_snowdp', 'f_zwt']
def means(tag, s):
    acc = {k: [] for k in V + ['t10', 't30', 'ch4']}
    for f in sorted(glob.glob(f'{B}sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        d = nc.Dataset(f)
        for k in V:
            acc[k].append(float(np.nanmean(np.ma.filled(d[k][:].astype(float), np.nan))))
        t = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)   # (time, lev, ...)
        t = t[:, 0, :]                                        # (time, soilsnow): 5 snow + 10 soil levels
        acc['t10'].append(float(np.nanmean(t[:, 5 + 3])))     # soil layer 4, node 0.12 m
        acc['t30'].append(float(np.nanmean(t[:, 5 + 5])))     # soil layer 6, node 0.37 m
        ft = f.replace('_hist_', '_hist_tracer_')
        try:
            dt = nc.Dataset(ft); acc['ch4'].append(float(np.nanmean(np.ma.filled(dt['f_methane_surf_flux_tot_active'][:].astype(float), np.nan))) * 16.043e3 * 86400)
        except Exception: pass
    return {k: np.nanmean(v) if v else np.nan for k, v in acc.items()}
for kind, a, b in (('main', 'v2_k4p_sasu5', 'v2_k4ph_sasu5'), ('val', 'v2_k4p_val', 'v2_k4ph_val')):
    sites = sorted(p.split('/')[-1] for p in glob.glob(f'{B}sites_{b}/*') if '_conf' not in p)
    print(f'== {kind}: site | ustar | H | LE | G | T~12cm | T~37cm | snow | zwt | CH4 mg/m2/d   (20 m -> tower height)')
    for s in sites:
        x, y = means(a, s), means(b, s)
        f = lambda k, fmt: f'{x[k]:{fmt}}->{y[k]:{fmt}}'
        print(f"{s:7s} {f('f_ustar','.2f')} {f('f_fsena','.0f')} {f('f_lfevpa','.0f')} {f('f_fgrnd','.1f')} {f('t10','.1f')} {f('t30','.1f')} {f('f_snowdp','.2f')} {f('f_zwt','.2f')} {f('ch4','.1f')}")
