"""Where the wetland-tile evapotranspiration comes from (Q-20): per tower, the
tower-height run's annual mean latent heat split into ground evaporation
(f_fevpg), canopy transpiration (f_etr) and canopy evaporation (f_fevpl - f_etr),
in W m-2, with the tower's observed latent heat and the mean LAI."""
import glob, re
import numpy as np, netCDF4 as nc
B = '/share/home/dq076/mode/Methane/cases/paper_v2/v260925u/sp_v2_c17c/'
LV = 2.501e6
print('site     LE_model | ground  transp  canopy_evap (W m-2, share) | LAI')
for kind, tag in (('main', 'v2_k4ph_sasu5'), ('val', 'v2_k4ph_val')):
    for p in sorted(glob.glob(f'{B}sites_{tag}/*')):
        s = p.split('/')[-1]
        if s == '_conf': continue
        acc = {k: [] for k in ('le', 'g', 't', 'l', 'lai')}
        for f in sorted(glob.glob(f'{p}/history/{s}_hist_2*.nc')):
            d = nc.Dataset(f)
            m = lambda v: float(np.nanmean(np.ma.filled(d[v][:].astype(float), np.nan)))
            acc['le'].append(m('f_lfevpa')); acc['g'].append(m('f_fevpg') * LV)
            acc['t'].append(m('f_etr') * LV); acc['l'].append((m('f_fevpl') - m('f_etr')) * LV); acc['lai'].append(m('f_lai'))
        if not acc['le']: continue
        a = {k: np.mean(v) for k, v in acc.items()}
        tot = a['g'] + a['t'] + a['l']
        print(f"{s:7s} {a['le']:6.0f} | {a['g']:5.0f} ({a['g']/tot:.2f}) {a['t']:5.0f} ({a['t']/tot:.2f}) {a['l']:5.0f} ({a['l']/tot:.2f}) | {a['lai']:.2f}")
