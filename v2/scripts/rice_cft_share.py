"""Rainfed (CFT 47 = PFT 61) and irrigated (CFT 48 = PFT 62) rice in the global
CFT surface data, and the irrigation method the runtime map gives each rice CFT
(the CH4 module treats a rice PFT as a flooded paddy only under method flood or
paddy). Q-15."""
import numpy as np, netCDF4 as nc
d = nc.Dataset('/share/home/dq076/data/CoLMrawdata_glwd/global_CFT_surface_data.nc')
print('CFT file vars:', [(v, d[v].shape) for v in d.variables][:8])
p = d['PCT_CFT']
lat = d['lat'][:]
w = np.cos(np.deg2rad(lat))[:, None]
pc = None
for name in ('PCT_CROP', 'pct_crop', 'PCTCROP'):
    if name in d.variables: pc = np.ma.filled(d[name][:].astype(float), 0); print('crop share var', name, pc.shape)
for k, lab in ((46, 'CFT47 rainfed rice'), (47, 'CFT48 irrigated rice')):
    x = np.ma.filled(p[:, :, k].astype(float), 0)
    print(lab, 'area-weighted sum of PCT_CFT (x cos lat):', round(float((x * w).sum()), 1))
m = nc.Dataset('/share/home/dq013/zhwei/colm/data/CoLMruntime/crop/surfdata_irrigation_method_96x144.nc')
print('irrigation method file vars:', [(v, m[v].shape) for v in m.variables])
cft = d['cft'][:]
print('cft coordinate at 47, 48 (1-based):', cft[46], cft[47])
for v in m.variables:
    if m[v].ndim >= 2:
        a = np.ma.filled(m[v][:], -1)
        print(v, m[v].dimensions, 'unique values:', np.unique(a)[:10])
# share of each rice CFT's area whose coarse cell carries each irrigation method
im = np.ma.filled(m['irrigation_method'][:].astype(float), -1)
mlat, mlon = m['lat'][:], m['lon'][:]
dlon = d['lon'][:]
ii = np.abs(lat[:, None] - mlat[None, :]).argmin(1)
jj = np.abs(((dlon[:, None] - mlon[None, :] + 180) % 360) - 180).argmin(1)
for k, cftid in ((46, 61), (47, 62)):
    x = np.ma.filled(p[:, :, k].astype(float), 0) * w
    meth = im[k][ii][:, jj]
    tot = x.sum()
    print(f'PFT {cftid}:', ', '.join(f'method {int(v)} {x[meth == v].sum() / tot:.3f}' for v in np.unique(meth)))
