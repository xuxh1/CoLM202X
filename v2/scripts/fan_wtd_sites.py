"""Fan et al. (2013) climatological water-table depth at the 23 main towers:
the 0.1-degree cell holding the tower and the 3x3 cells around it, annual
mean and the monthly minimum and maximum [m below surface], with the tower's
observed median (K4 summary). Water-table options for the wetland tile (Q-18)."""
import csv, re
import numpy as np, netCDF4 as nc
d = nc.Dataset('/share/home/dq013/zhwei/colm/data/CoLMruntime/wtd_fan2020_3600x1800.nc')
lat, lon = d['lat'][:], d['lon'][:]
S = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Sitedata/'
summ = open('/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/results/v2_k4_sasu5_summary.txt').read()
obs = {m.group(1): float(m.group(2)) for m in re.finditer(r'^(\w\w-\w\w\w)\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+\d', summ, re.M)}
cls = {m.group(1): m.group(2) for m in re.finditer(r'^(\w\w-\w\w\w)\s+(\w+)\s+-?\d', summ, re.M)}
import glob
print('site    class        obs_med  fan_cell_mean fan_cell_min fan_cell_max  fan_3x3_min_of_means')
for f in sorted(glob.glob(S + '*_site.nc')):
    sid = f.split('/')[-1][:6]
    if sid not in cls: continue
    s = nc.Dataset(f); la, lo = float(s['latitude'][...]), float(s['longitude'][...])
    i = int(np.argmin(abs(lat - la))); j = int(np.argmin(abs(((lon - lo + 180) % 360) - 180)))
    w = np.ma.filled(d['wtd'][:, i, j].astype(float), np.nan)
    blk = np.ma.filled(d['wtd'][:, max(i-1,0):i+2, max(j-1,0):j+2].astype(float), np.nan).mean(axis=0)
    print(f'{sid} {cls[sid]:11s} {obs.get(sid, np.nan):7.2f} {np.nanmean(w):12.2f} {np.nanmin(w):12.2f} {np.nanmax(w):12.2f} {np.nanmin(blk):14.2f}')
