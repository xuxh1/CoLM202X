"""Winter screen over every tower with a frozen season (Q-27): whether the
model soil is too cold in the three coldest months and whether a missing
snowpack goes with it. Per tower, over the years of the tower CH4 record:
forcing precipitation of the months with mean air temperature below 0 C
(mm per cold season) and the share of it that falls at air temperature above
0 C and below -2 C (scheme II phase split turns the first into rain); model
and tower snow depth and soil temperature at the tower probe nearest 10 cm
(depths from the FLUXNET-CH4 META table) in the three coldest months; and
November-to-March CH4: share of half-hours with a measured flux, mean of the
measured half-hours, gap-filled mean, model mean.
Usage: winter_screen.py <case under cases/> <tag> [<tag> ...]"""
import glob
import os
import re
import sys
from collections import defaultdict

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cold_obs as CO                                                # noqa: E402

M = CO.M
FORC = f'{M}/data/FLUXNET-CH4/Forcing'
NDJFM = (11, 12, 1, 2, 3)


def forcing(sid, yrs):
    p = glob.glob(f'{FORC}/{sid}_*_Met.nc')[0]
    with nc.Dataset(p) as d:
        t = nc.num2date(d['time'][:], d['time'].units)
        ta = np.ma.filled(d['Tair'][:].astype(float), np.nan).ravel() - 273.15
        pr = np.ma.filled(d['Precip'][:].astype(float), np.nan).ravel()
        dt = float(d['time'][1] - d['time'][0])
    ym = np.array([(x.year, x.month) for x in t])
    keep = np.isin(ym[:, 0], yrs)
    tm = {m: np.nanmean(ta[keep & (ym[:, 1] == m)]) for m in range(1, 13)}
    cold = [m for m in range(1, 13) if tm[m] < 0.]
    coldest = sorted(range(1, 13), key=lambda m: tm[m])[:3]
    sel = keep & np.isin(ym[:, 1], cold)
    tot = np.nansum(pr[sel]) * dt
    nyr = len(set(ym[keep, 0]))
    warm = np.nansum(pr[sel & (ta > 0.)]) * dt / tot if tot > 0 else np.nan
    vcold = np.nansum(pr[sel & (ta < -2.)]) * dt / tot if tot > 0 else np.nan
    return tot / max(nyr, 1), warm, vcold, coldest, tm


def main():
    case, tags = sys.argv[1], sys.argv[2:]
    dep = CO.probes()
    print('site    tag       | cold-P mm  >0C  <-2C | snow mod/obs m | Tsoil@z  z(m)  obs    mod   | NDJFM meas.share  meas  gapf   mod')
    for tag in tags:
        for sdir in sorted(glob.glob(f'{M}/cases/{case}/sites_{tag}/[A-Z][A-Z]-*')):
            sid = os.path.basename(sdir)
            try:
                t, v = CO.tower(sid)
            except IndexError:
                continue
            ym = np.array([(x.year, x.month) for x in t])
            yrs = sorted({y for y, m in ym[np.isfinite(v['FCH4_f'])]})
            try:
                cp, warm, vcold, coldest, tm = forcing(sid, yrs)
            except (IndexError, KeyError):
                continue
            if tm[coldest[0]] >= 0.:
                continue                                          # no frozen season
            # tower winter values
            wsel = np.isin(ym[:, 1], coldest)
            osnow = np.nanmean(v['SnowDepth'][wsel]) if v['SnowDepth'] is not None and np.isfinite(v['SnowDepth'][wsel]).any() else np.nan
            pz = {k: z for k, z in dep.get(sid, {}).items() if k <= 3 and z > 0}
            k10 = min(pz, key=lambda k: abs(pz[k] - 0.1)) if pz else None
            ots = np.nanmean(v[f'TS_{k10}'][wsel]) if k10 and np.isfinite(v[f'TS_{k10}'][wsel]).any() else np.nan
            csel = np.isin(ym[:, 1], NDJFM)
            share = np.isfinite(v['FCH4'][csel]).mean()
            meas = np.nanmean(v['FCH4'][csel]) * CO.TO_MG if np.isfinite(v['FCH4'][csel]).any() else np.nan
            gapf = np.nanmean(v['FCH4_f'][csel]) * CO.TO_MG
            # model
            ms, mt, me = [], [], []
            for y in yrs:
                try:
                    d = nc.Dataset(f'{sdir}/history/{sid}_hist_{y}.nc')
                    tr = nc.Dataset(f'{sdir}/history/{sid}_hist_tracer_{y}.nc')
                except OSError:
                    continue
                with d, tr:
                    T = np.ma.filled(d['f_t_soisno'][:].astype(float), np.nan)[:, 0, 5:15] - 273.15
                    if T.shape[0] < 12:
                        continue
                    sd = np.ma.filled(d['f_snowdp'][:].astype(float), np.nan).reshape(12, -1)[:, 0]
                    e = np.ma.filled(tr['f_methane_surf_flux_tot_active'][:].astype(float), np.nan).reshape(12, -1)[:, 0] * 16.043e3 * 86400
                    for m in coldest:
                        ms.append(sd[m - 1])
                        if k10:
                            mt.append(np.interp(pz[k10], CO.ZSOI, T[m - 1]))
                    for m in NDJFM:
                        me.append(e[m - 1])
            f = lambda x, w=6, p=1: f'{x:{w}.{p}f}' if np.isfinite(x) else ' ' * (w - 1) + '.'
            print(f'{sid:7s} {tag[-5:]:9s} | {f(cp, 7, 0)} {f(warm, 5, 2)} {f(vcold, 5, 2)} | '
                  f'{f(np.mean(ms) if ms else np.nan, 5, 2)} {f(osnow, 5, 2)} | '
                  f'{f(pz[k10], 5, 3) if k10 else "    ."} {f(ots)} {f(np.mean(mt) if mt else np.nan)} | '
                  f'{f(share, 5, 2)} {f(meas)} {f(gapf)} {f(np.mean(me) if me else np.nan)}')


if __name__ == '__main__':
    main()
