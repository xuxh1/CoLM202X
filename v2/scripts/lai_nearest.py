"""Nearest pixels with a valid LAI around a tower in the 15-arcsecond raw data,
by IGBP class (US-WPT, US-ORv read LAI 0): for each class, the distance of the
nearest pixel whose peak monthly LAI exceeds 0.5 and the mean peak of the
nearest four. Usage: lai_nearest.py <lat> <lon> <year> [half-width pixels]"""
import sys
import numpy as np, netCDF4 as nc
R = '/share/home/dq013/zhwei/colm/data/CoLMrawdata'
lat, lon, year = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
hw = int(sys.argv[4]) if len(sys.argv) > 4 else 40
la0, lo0 = int(np.ceil(lat / 5) * 5), int(np.floor(lon / 5) * 5)
with nc.Dataset(f'{R}/plant_15s/RG_{la0}_{lo0}_{la0 - 5}_{lo0 + 5}.MOD{year}.nc') as d:
    la, lo = d['lat'][:], d['lon'][:]
    i, j = int(np.argmin(abs(la - lat))), int(np.argmin(abs(lo - lon)))
    sl = (slice(max(i - hw, 0), i + hw + 1), slice(max(j - hw, 0), j + hw + 1))
    lc = np.asarray(d['LC'][sl])
    pk = np.ma.filled(d['MONTHLY_LC_LAI'][(slice(None),) + sl].astype(float), 0).max(0)
    y = (la[sl[0]] - lat)[:, None] * 111.2
    x = (lo[sl[1]] - lon)[None, :] * 111.2 * np.cos(np.radians(lat))
dist = np.hypot(y, x)
print(f'tower pixel class {lc[i - sl[0].start, j - sl[1].start]}, window +-{hw} px')
for c in np.unique(lc):
    m = (lc == c) & (pk > 0.5)
    if not m.any():
        print(f'  class {c:2d}: {int((lc == c).sum()):5d} px, none with LAI')
        continue
    o = np.argsort(dist[m])[:4]
    print(f'  class {c:2d}: {int((lc == c).sum()):5d} px, {int(m.sum()):5d} with LAI; nearest {dist[m][o[0]]:.2f} km, '
          f'nearest-4 mean peak {pk[m][o].mean():.1f}')
