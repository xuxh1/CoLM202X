"""Land class and peak LAI of the tower pixel year by year in the 15-arcsecond
raw data (zero-LAI tower years). Usage: lai_pixel_years.py <lat> <lon> <y0> <y1>"""
import sys
import numpy as np, netCDF4 as nc
R = '/share/home/dq013/zhwei/colm/data/CoLMrawdata/plant_15s'
lat, lon, y0, y1 = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
la0, lo0 = int(np.ceil(lat / 5) * 5), int(np.floor(lon / 5) * 5)
for y in range(y0, y1 + 1):
    with nc.Dataset(f'{R}/RG_{la0}_{lo0}_{la0 - 5}_{lo0 + 5}.MOD{y}.nc') as d:
        i, j = int(np.argmin(abs(d['lat'][:] - lat))), int(np.argmin(abs(d['lon'][:] - lon)))
        lai = np.ma.filled(d['MONTHLY_LC_LAI'][:, i, j].astype(float), 0.)
        print(y, 'class', int(d['LC'][i, j]), 'peak LAI', round(float(lai.max()), 2), 'PCT_WETLAND', float(d['PCT_WETLAND'][i, j]))
