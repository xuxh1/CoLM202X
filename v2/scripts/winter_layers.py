#!/usr/bin/env python
"""Where does a tower's winter CH4 come from? Monthly climatology per soil
layer: temperature (C), ice share of the layer's water, saturated-subcolumn
production (mg CH4 m-2 d-1 from that layer) and soil plus litter carbon
(kg C m-2 in that layer). Usage: winter_layers.py <case under cases/> <tag> <site> [months]"""
import glob
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd

R = '/share/home/dq076/mode/Methane'
MG_D = 16.04e3 * 86400.0                                            # mol s-1 -> mg d-1
POOLS = ('f_litr1c_vr', 'f_litr2c_vr', 'f_litr3c_vr', 'f_soil1c_vr', 'f_soil2c_vr', 'f_soil3c_vr')


def layers(nl):
    """CoLM node depths and thicknesses (main/MOD_Vars_Global.F90:138-147); the
    history 'soil' coordinate is only the layer index."""
    z = 0.025 * (np.exp(0.5 * (np.arange(1, nl + 1) - 0.5)) - 1.0)
    dz = np.empty(nl)
    dz[0] = 0.5 * (z[0] + z[1])
    dz[-1] = z[-1] - z[-2]
    dz[1:-1] = 0.5 * (z[2:] - z[:-2])
    return z, dz


def main():
    case, tag, s = sys.argv[1:4]
    months = [int(m) for m in sys.argv[4].split(',')] if len(sys.argv) > 4 else [10, 11, 12, 1, 2, 3, 4, 5, 7]
    acc = {}
    for f in sorted(glob.glob(f'{R}/cases/{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        with nc.Dataset(f) as h, nc.Dataset(f.replace('_hist_', '_hist_tracer_')) as e:
            if h['f_t_soisno'].shape[0] != 12:
                continue
            g = lambda d, v: np.ma.filled(d[v][:], np.nan).astype(float)[:, 0]
            nl = h.dimensions['soil'].size
            z, dz = layers(nl)
            T = g(h, 'f_t_soisno')[:, -nl:] - 273.15
            wl, wi = g(h, 'f_wliq_soisno')[:, -nl:], g(h, 'f_wice_soisno')[:, -nl:]
            P = g(e, 'f_methane_prod_depth_sat')[:, :nl] * dz * MG_D
            C = sum(g(h, v)[:, :nl] for v in POOLS) * dz / 1e3
            for k, a in (('T', T), ('ice', wi / np.maximum(wl + wi, 1e-9)), ('P', P), ('C', C)):
                acc.setdefault(k, []).append(a)
    print(f'== {s} ({tag}); layer node depths (m): ' + ' '.join(f'{v:.2f}' for v in z))
    for k, lab in (('T', 'temperature C'), ('ice', 'ice share'), ('P', 'sat production mg m-2 d-1'),
                   ('C', 'soil+litter C kg m-2')):
        a = np.nanmean(np.stack(acc[k]), axis=0)                    # 12 x nl
        df = pd.DataFrame(a[[m - 1 for m in months]], index=[f'm{m}' for m in months],
                          columns=[f'L{j + 1}' for j in range(nl)])
        if k == 'P':
            df['column'] = df.sum(axis=1)
        print(f'-- {lab}')
        print(df.round(3 if k == 'ice' else 2).to_string())


if __name__ == '__main__':
    main()
