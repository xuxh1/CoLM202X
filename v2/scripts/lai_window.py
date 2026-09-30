"""LAI and land type around a tower in the 15-arcsecond raw data the site runs
extract from (US-WPT, US-ORv read water pixels): a window of pixels centred on
the tower, each with its IGBP class and the year's peak and mean monthly LAI.
Usage: lai_window.py <lat> <lon> <year> [half-width in pixels]"""
import glob, sys
import numpy as np, netCDF4 as nc
R = '/share/home/dq013/zhwei/colm/data/CoLMrawdata'
lat, lon, year = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
hw = int(sys.argv[4]) if len(sys.argv) > 4 else 4
f = f'{R}/plant_15s/RG_{int(np.ceil(lat/5)*5)}_{int(np.floor(lon/5)*5)}_{int(np.ceil(lat/5)*5)-5}_{int(np.floor(lon/5)*5)+5}.MOD{year}.nc'
with nc.Dataset(f) as d:
    print(f.split('/')[-1], list(d.variables))
    la, lo = d['lat'][:], d['lon'][:]
    i, j = int(np.argmin(abs(la - lat))), int(np.argmin(abs(lo - lon)))
    lt = d['LC'][i - hw:i + hw + 1, j - hw:j + hw + 1]
    pw = np.ma.filled(d['PCT_WATER'][i - hw:i + hw + 1, j - hw:j + hw + 1].astype(float), np.nan)
    v = 'MONTHLY_LC_LAI'
    lai = np.ma.filled(d[v][..., i - hw:i + hw + 1, j - hw:j + hw + 1].astype(float), np.nan)
    print('LAI var', v, d[v].dimensions, lai.shape)
    lai = lai.reshape(-1, lai.shape[-2], lai.shape[-1])
print('rows north to south as stored; centre = tower; cell = IGBP:peak/mean LAI:water %')
for r in range(2 * hw + 1):
    print(' '.join(f'{int(lt[r, c]):2d}:{np.nanmax(lai[:, r, c]):3.1f}/{np.nanmean(lai[:, r, c]):3.1f}:{pw[r, c]:3.0f}' for c in range(2 * hw + 1)))
