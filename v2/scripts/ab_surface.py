"""Surface exchange and soil state of two site tags against each other and the
towers (V-36): per tower the mean over shared months of friction velocity,
sensible and latent heat and net radiation (tower Ustar, Qh_f, Qle_f, Rnet_f
against f_ustar, f_fsena, f_lfevpa, f_rnet) and the evaporative fraction
LE/(H + LE) over the months all three have, then the annual means of ground evaporation, snow depth, soil
temperature near 0.12 m, water-table depth and CH4 flux of both tags; last the
median ratio model/obs over towers for A and B.
Usage: ab_surface.py <case A> <tag A> <case B> <tag B>  (cases under cases/)"""
import glob
import re
import sys
from collections import defaultdict
import numpy as np
import netCDF4 as nc

C = '/share/home/dq076/mode/Methane/cases/'
OBS = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Observation/'
PAIRS = {'ustar': ('Ustar', 'f_ustar'), 'H': ('Qh_f', 'f_fsena'), 'LE': ('Qle_f', 'f_lfevpa'),
         'Rn': ('Rnet_f', 'f_rnet')}
LV = 2.501e6


def obs_monthly(s):
    out = {k: defaultdict(list) for k in PAIRS}
    for p in glob.glob(f'{OBS}{s}_*_Flux.nc'):
        with nc.Dataset(p) as d:
            t = nc.num2date(d['time'][:], d['time'].units)
            keys = [(x.year, x.month) for x in t]
            for k, (ov, _) in PAIRS.items():
                if ov not in d.variables:
                    continue
                v = np.ma.filled(d[ov][:].astype(float), np.nan).ravel()
                for kk, vv in zip(keys, v):
                    out[k][kk].append(vv)
    return {k: {kk: np.nanmean(vv) for kk, vv in m.items() if np.isfinite(vv).sum() >= 240} for k, m in out.items()}


def mod(case, tag, s):
    """(monthly surface fluxes by key, annual means of the state variables)"""
    mon = {k: {} for k in PAIRS}
    acc = defaultdict(list)
    for f in sorted(glob.glob(f'{C}{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        y = int(re.search(r'_(\d{4})\.nc$', f).group(1))
        with nc.Dataset(f) as d:
            g = lambda v: np.ma.filled(d[v][:].astype(float), np.nan).reshape(d[v].shape[0], -1)[:, 0]
            for k, (_, mv) in PAIRS.items():
                v = g(mv)
                for i in range(min(12, v.size)):
                    mon[k][(y, i + 1)] = v[i]
            acc['Eg'].append(np.nanmean(g('f_fevpg')) * LV)
            acc['snow'].append(np.nanmean(g('f_snowdp')))
            acc['zwt'].append(np.nanmean(g('f_zwt')))
            t = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)[:, 0, :]
            acc['t12'].append(np.nanmean(t[:, 5 + 3]))          # soil layer 4, node 0.12 m
        ft = f.replace('_hist_', '_hist_tracer_')
        try:
            with nc.Dataset(ft) as d:
                acc['ch4'].append(np.nanmean(np.ma.filled(d['f_methane_surf_flux_tot_active'][:].astype(float), np.nan))
                                  * 16.043e3 * 86400)
        except (OSError, KeyError):
            pass
    return mon, {k: np.nanmean(v) for k, v in acc.items()}


def main():
    ca, ta, cb, tb = sys.argv[1:5]
    sites = sorted(p.split('/')[-1] for p in glob.glob(f'{C}{cb}/sites_{tb}/*') if '_conf' not in p)
    ratio = {k: {'a': [], 'b': []} for k in list(PAIRS) + ['EF']}
    print(f'# A = {ta} ({ca}); B = {tb} ({cb})')
    print('site    | ustar obs/A/B | H obs/A/B | LE obs/A/B | Rn obs/A/B (W m-2) | EF obs/A/B | Eg A/B | snow m A/B | T12 K A/B | zwt m A/B | CH4 mg m-2 d-1 A/B')
    for s in sites:
        o = obs_monthly(s)
        ma, sa = mod(ca, ta, s)
        mb, sb = mod(cb, tb, s)
        cells = []
        for k in PAIRS:
            ks = [kk for kk in o[k] if kk in ma[k] and kk in mb[k]]
            if len(ks) < 6:
                cells.append('   n/a')
                continue
            O = np.mean([o[k][kk] for kk in ks])
            A = np.mean([ma[k][kk] for kk in ks])
            B = np.mean([mb[k][kk] for kk in ks])
            f = '.2f' if k == 'ustar' else '.0f'
            cells.append(f'{O:{f}}/{A:{f}}/{B:{f}}')
            if abs(O) > 1e-3:
                ratio[k]['a'].append(A / O)
                ratio[k]['b'].append(B / O)
        ks = [kk for kk in o['H'] if kk in o['LE'] and kk in ma['H'] and kk in mb['H']]
        if len(ks) >= 6:
            ef = []
            for src in (o, ma, mb):
                h = np.mean([src['H'][kk] for kk in ks]); le = np.mean([src['LE'][kk] for kk in ks])
                ef.append(le / (h + le) if h + le > 1. else np.nan)
            cells.append('/'.join(f'{x:.2f}' for x in ef))
            if np.isfinite(ef).all():
                ratio['EF']['a'].append(ef[1] / ef[0]); ratio['EF']['b'].append(ef[2] / ef[0])
        else:
            cells.append('   n/a')
        st = (f"{sa.get('Eg', np.nan):3.0f}/{sb.get('Eg', np.nan):3.0f} | {sa.get('snow', np.nan):.2f}/{sb.get('snow', np.nan):.2f} | "
              f"{sa.get('t12', np.nan):.1f}/{sb.get('t12', np.nan):.1f} | {sa.get('zwt', np.nan):.2f}/{sb.get('zwt', np.nan):.2f} | "
              f"{sa.get('ch4', np.nan):5.1f}/{sb.get('ch4', np.nan):5.1f}")
        print(f'{s:7s} | ' + ' | '.join(cells) + ' | ' + st)
    print('# median model/obs over towers: A | B')
    for k in ratio:
        print(f'{k:6s} {np.median(ratio[k]["a"]):.2f} | {np.median(ratio[k]["b"]):.2f}  (n={len(ratio[k]["a"])})')


if __name__ == '__main__':
    main()
