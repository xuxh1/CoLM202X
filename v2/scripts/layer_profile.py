"""Where in the column the summer CH4 comes from at the towers (main-set
diagnosis): June-August (December-February south of the equator) climatology
of the layer shares of fresh carbon (litter plus soil1 and soil2), potential
decomposition CO2, CH4 production and oxidation, with root water uptake as
the root profile, the water-table depth and the inundated fraction. Layer
interfaces are CoLM's default soil levels.
Usage: layer_profile.py <case under cases/> <tag> <site> [site ...]"""
import glob
import sys
import numpy as np
import netCDF4 as nc

C = '/share/home/dq076/mode/Methane/cases/'
ZI = [0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289, 1.3828, 2.2961, 3.8019]
NL = 8


def clim(case, tag, s, months):
    acc = {}
    for f in sorted(glob.glob(f'{C}{case}/sites_{tag}/{s}/history/{s}_hist_2*.nc')):
        ft = f.replace('_hist_', '_hist_tracer_')
        with nc.Dataset(f) as d, nc.Dataset(ft) as t:
            lat = float(np.ravel(d['lat'][:])[0]) if 'lat' in d.variables else 60.
            mo = [m - 1 for m in (months if lat >= 0 else [12, 1, 2])]
            g = lambda ds, v: np.ma.filled(ds[v][:].astype(float), np.nan)[mo].reshape(len(mo), -1)
            c = sum(g(d, v)[:, :NL] for v in ('f_litr1c_vr', 'f_litr2c_vr', 'f_litr3c_vr', 'f_soil1c_vr', 'f_soil2c_vr'))
            for k, a in (('C', c), ('root', g(d, 'f_rootr')[:, :NL]), ('decomp', g(t, 'f_co2_decomp_depth')[:, :NL]),
                         ('prod', g(t, 'f_methane_prod_depth')[:, :NL]), ('oxid', g(t, 'f_methane_oxid_depth')[:, :NL])):
                acc.setdefault(k, []).append(np.nanmean(a, 0))
            acc.setdefault('zwt', []).append(np.nanmean(g(d, 'f_zwt')))
            try:
                acc.setdefault('finund', []).append(np.nanmean(g(t, 'f_methane_finundated')))
            except (IndexError, KeyError):
                pass
    return {k: np.nanmean(np.array(v), 0) for k, v in acc.items()}


def main():
    case, tag = sys.argv[1:3]
    dz = np.diff([0.] + ZI)[:NL]
    print(f'# {tag}: summer layer shares (%); layer bottoms (m) ' + ' '.join(f'{z:.2f}' for z in ZI[:NL]))
    for s in sys.argv[3:]:
        a = clim(case, tag, s, [6, 7, 8])
        print(f'{s}: zwt {a["zwt"]:.2f} m, inundated fraction {a.get("finund", np.nan):.2f}')
        for k in ('C', 'root', 'decomp', 'prod', 'oxid'):
            v = a[k] * (dz if k in ('decomp', 'prod', 'oxid', 'C') else 1.)
            sh = 100 * v / np.nansum(v) if np.nansum(v) > 0 else v * np.nan
            print(f'   {k:6s} ' + ' '.join(f'{x:4.0f}' for x in sh))
        p, o = np.nansum(a['prod'] * dz), np.nansum(a['oxid'] * dz)
        print(f'   oxid/prod {o / p:.2f}' if p > 0 else '   no production')


if __name__ == '__main__':
    main()
