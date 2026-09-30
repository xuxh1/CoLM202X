"""Surface exchange against the towers (Q-20): monthly friction velocity, sensible
and latent heat from the tower (Ustar, Qh_f, Qle_f) against the 20 m run
(v2_k4p) and the tower-height run (v2_k4ph); per tower the mean over the months
both have, and the median ratio model/obs over towers."""
import glob, re
from collections import defaultdict
import numpy as np, netCDF4 as nc
B = '/share/home/dq076/mode/Methane/cases/paper_v2/v260925u/sp_v2_c17c/'
OBS = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Observation/'
PAIRS = {'ustar': ('Ustar', 'f_ustar'), 'H': ('Qh_f', 'f_fsena'), 'LE': ('Qle_f', 'f_lfevpa')}
def obs_monthly(s):
    out = {k: defaultdict(list) for k in PAIRS}
    for p in glob.glob(f'{OBS}{s}_*_Flux.nc'):
        d = nc.Dataset(p)
        t = nc.num2date(d['time'][:], d['time'].units)
        keys = [(x.year, x.month) for x in t]
        for k, (ov, _) in PAIRS.items():
            if ov not in d.variables: continue
            v = np.ma.filled(d[ov][:].astype(float), np.nan).ravel()
            for kk, vv in zip(keys, v): out[k][kk].append(vv)
    return {k: {kk: np.nanmean(vv) for kk, vv in m.items() if np.isfinite(vv).sum() >= 240} for k, m in out.items()}
def mod_monthly(tag, s):
    out = {k: {} for k in PAIRS}
    for f in sorted(glob.glob(f'{B}sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        y = int(re.search(r'_(\d{4})\.nc$', f).group(1)); d = nc.Dataset(f)
        for k, (_, mv) in PAIRS.items():
            v = np.ma.filled(d[mv][:].astype(float), np.nan).ravel()
            for i in range(min(12, v.size)): out[k][(y, i + 1)] = v[i]
    return out
res = {k: {'a': [], 'b': []} for k in PAIRS}
for kind, a, b in (('main', 'v2_k4p_sasu5', 'v2_k4ph_sasu5'), ('val', 'v2_k4p_val', 'v2_k4ph_val')):
    sites = sorted(p.split('/')[-1] for p in glob.glob(f'{B}sites_{b}/*') if '_conf' not in p)
    print(f'== {kind}: site | ustar obs / 20m / tower | H obs / 20m / tower | LE obs / 20m / tower  (W m-2, m s-1; mean over shared months)')
    for s in sites:
        o, ma, mb = obs_monthly(s), mod_monthly(a, s), mod_monthly(b, s)
        cells = []
        for k in PAIRS:
            ks = [kk for kk in o[k] if kk in ma[k] and kk in mb[k]]
            if len(ks) < 6: cells.append('   n/a'); continue
            O = np.mean([o[k][kk] for kk in ks]); A = np.mean([ma[k][kk] for kk in ks]); Bv = np.mean([mb[k][kk] for kk in ks])
            fmt = '.2f' if k == 'ustar' else '.0f'
            cells.append(f'{O:{fmt}}/{A:{fmt}}/{Bv:{fmt}}')
            if abs(O) > 1e-3:
                res[k]['a'].append(A / O); res[k]['b'].append(Bv / O)
        print(f'{s:7s} ' + ' | '.join(cells))
print('== median model/obs over towers: 20 m | tower height')
for k in PAIRS:
    print(f'{k:6s} {np.median(res[k]["a"]):.2f} | {np.median(res[k]["b"]):.2f}  (n={len(res[k]["a"])})')
