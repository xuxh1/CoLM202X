"""scripts/sites/score_series.py with the V2 observation gate (ledger D-17).

The shared scorer counts finite half-hours when it gates a month (>= 10 of
them, about five hours), although its docstring asks for 10 valid days, so
months with one or two days of data enter the KGE terms (Q-16: 38 site-months
at 15 of the 23 towers). Here a day is valid when at least half of its time
steps are finite, its value is the mean of those steps, and a month needs 10
valid days; its value is the mean of its daily values. Segments the Q-16
audit found unfit as a target are left out:
  DE-Hte 2011-04  start-up spikes up to 31921 nmol m-2 s-1 (monthly 785 against 78-168)
  US-Sne 2016     pasture before the November 2016 rewetting
Everything else, output files and command line included, is score_series'.
Usage: score_series_v2.py <tag> <case dir> [--sitelist csv]"""
import glob
import sys
from collections import defaultdict

import numpy as np
import netCDF4 as nc

sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/sites')
import score_series as S                                            # noqa: E402

EXCLUDE = {'DE-Hte': lambda y, m: (y, m) == (2011, 4),
           'US-Sne': lambda y, m: y == 2016}
MIN_DAYS = 10


def obs_monthly(sid):
    """(year, month) -> mean of valid daily means, >= MIN_DAYS valid days."""
    days = defaultdict(list)
    for p in glob.glob(f'{S.OBS}/{sid}_*_Flux.nc'):
        with nc.Dataset(p) as d:
            if 'FCH4_f' not in d.variables:
                continue
            t = nc.num2date(d.variables['time'][:], d.variables['time'].units)
            v = np.ma.filled(d.variables['FCH4_f'][:], np.nan).ravel() * S.TO_MGM2D
            for ti, vi in zip(t, v):
                days[(ti.year, ti.month, ti.day)].append(vi)
    months = defaultdict(list)
    for (y, m, _), vs in days.items():
        a = np.asarray(vs, dtype=float)
        if np.isfinite(a).sum() >= 0.5 * a.size:
            months[(y, m)].append(float(np.nanmean(a)))
    skip = EXCLUDE.get(sid, lambda y, m: False)
    return {k: float(np.mean(v)) for k, v in months.items()
            if len(v) >= MIN_DAYS and not skip(*k)}


S.obs_monthly = obs_monthly

if __name__ == '__main__':
    S.main()
