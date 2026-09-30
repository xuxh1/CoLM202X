"""Where a moss / peat surface resistance would act (Q-20, C-26 design): per
tower of the tower-height run, the ground evaporation of snow-free months
split by the modelled water level of that month (ponded above 1 mm, table
0-0.10 m, 0.10-0.30 m, deeper), with the resistance r_w(z) = 100 s m-1 up to
0.10 m and min(1000, 100 exp((z - 0.10)/0.155)) below, and the implied
ground evaporation factor ra/(ra + r_w) for a typical ra of 50 s m-1."""
import glob
import numpy as np
import netCDF4 as nc
B = '/share/home/dq076/mode/Methane/cases/paper_v2/v260925u/sp_v2_c17c/'
LV = 2.501e6
RA = 50.


def rw(z):
    return np.where(z <= 0.10, 100., np.minimum(1000., 100. * np.exp((z - 0.10) / 0.155)))


print('site     Eg W/m2 | share of snow-free Eg: pond  0-0.1  0.1-0.3  >0.3 | Eg factor')
for kind, tag in (('main', 'v2_k4ph_sasu5'), ('val', 'v2_k4ph_val')):
    print(f'== {kind} {tag}')
    for p in sorted(glob.glob(f'{B}sites_{tag}/*')):
        s = p.split('/')[-1]
        if s == '_conf':
            continue
        eg, zw, pd, sn = [], [], [], []
        for f in sorted(glob.glob(f'{p}/history/{s}_hist_2*.nc')):
            with nc.Dataset(f) as d:
                g = lambda v: np.ma.filled(d[v][:].astype(float), np.nan).reshape(d[v].shape[0], -1)[:, 0]
                eg.append(g('f_fevpg') * LV); zw.append(g('f_zwt')); pd.append(g('f_wdsrf')); sn.append(g('f_fsno'))
        if not eg:
            continue
        eg, zw, pd, sn = (np.concatenate(x) for x in (eg, zw, pd, sn))
        e = eg * (1. - sn)                        # snow-free part of ground evaporation
        pond = pd > 1.
        bins = [pond, ~pond & (zw <= 0.10), ~pond & (zw > 0.10) & (zw <= 0.30), ~pond & (zw > 0.30)]
        tot = np.nansum(np.clip(e, 0, None))
        sh = [np.nansum(np.clip(e[b], 0, None)) / tot if tot > 0 else np.nan for b in bins]
        fac = np.where(pond, 1., RA / (RA + rw(zw)))
        f = np.nansum(np.clip(e, 0, None) * fac) / tot if tot > 0 else np.nan
        print(f'{s:7s} {np.nanmean(eg):6.1f} | {sh[0]:5.2f} {sh[1]:6.2f} {sh[2]:7.2f} {sh[3]:5.2f} | {f:5.2f}')
