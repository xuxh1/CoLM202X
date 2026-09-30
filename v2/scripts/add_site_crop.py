#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Give three drained-peat crop towers of the Sacramento-San Joaquin Delta the
crop actually grown under the tower (user 2026-09-30, Q-80).

Single-point mode builds a CROPLAND site from the coarse global CFT map unless
the site file carries croptyp and pctcrop and the run sets USE_SITE_pctcrop
(mksrfdata/MOD_SingleSrfdata.F90), so US-Bi1, US-Bi2 and US-Tw3 each got the
crops of their 2-degree cell, 10-14 % of it irrigated rice, and the model put
their CH4 in May-September from that rice (v2/results/diag_drained_260930.txt).
Each tower gets one crop patch, as scripts/sites/add_site_croptyp.py did for
the seven rice towers (CFT index, PFT = CFT + 14; main/MOD_Const_PFT.F90):
  US-Bi1  CFT 2  c3_irrigated              perennial alfalfa, flood irrigated about
          once a month May-September (Anthony and Silver 2023, Table 1, p. 465)
  US-Tw3  CFT 2  c3_irrigated              perennial alfalfa since 2010 (Hemes et al.
          2018, Sect. 2.1); irrigation not stated for this field, taken as the
          Delta alfalfa practice documented at US-Bi1 (inferred)
  US-Bi2  CFT 3  temperate_corn (rainfed)  corn, 118 kg N/ha/yr, winter flooding
          November-March (Anthony and Silver 2023, p. 464); summer irrigation not
          stated, so rainfed; the winter flooding is not represented by the crop
The Sitedata files are rewritten in place (one copy of the data only); md5 sums
before the change go to data/_tmp/分析_260930_塔气温与作物/. Without --apply only prints.
Usage: add_site_crop.py [--apply]"""
import datetime as dt
import glob
import hashlib
import os
import sys

import netCDF4 as nc
import numpy as np

SITEDATA = '/share/home/dq076/mode/Methane/data/FLUXNET-CH4/Sitedata'
LOG = '/share/home/dq076/mode/Methane/data/_tmp/分析_260930_塔气温与作物'
CROP = {'US-Bi1': (2, 'irrigated C3 crop (alfalfa, flood irrigated monthly May-Sep)'),
        'US-Tw3': (2, 'irrigated C3 crop (alfalfa, Delta practice as at US-Bi1)'),
        'US-Bi2': (3, 'rainfed temperate corn (irrigation not documented)')}


def main():
    apply = '--apply' in sys.argv[1:]
    md5 = []
    for sid, (cft, what) in CROP.items():
        ps = glob.glob(f'{SITEDATA}/{sid}_*_site.nc')
        assert len(ps) == 1, ps
        p = ps[0]
        with nc.Dataset(p) as d:
            igbp = int(d['IGBP_classification'][...])
            has = [v for v in ('croptyp', 'pctcrop') if v in d.variables]
        print(f'{sid}: {os.path.basename(p)}  IGBP {igbp}  existing {has or "none"}  -> CFT {cft} ({what})')
        assert igbp == 12, f'{sid} is not CROPLAND'
        if not apply:
            continue
        md5.append(f'{hashlib.md5(open(p, "rb").read()).hexdigest()}  {p}')
        with nc.Dataset(p, 'a') as d:
            if 'crop' not in d.dimensions:
                d.createDimension('crop', 1)
            v = d['croptyp'] if 'croptyp' in d.variables else d.createVariable('croptyp', 'i4', ('crop',))
            v[:] = np.array([cft], dtype='i4')
            v.long_name = 'crop functional type of the tower (CFT index, PFT = CFT + 14)'
            v.comment = f'{what}; read when USE_SITE_pctcrop is on'
            w = d['pctcrop'] if 'pctcrop' in d.variables else d.createVariable('pctcrop', 'f8', ('crop',))
            w[:] = np.array([1.0])
            w.long_name = 'share of the crop functional type in the tower footprint'
            w.units = '1'
            stamp = dt.datetime.now().strftime('%Y-%m-%d %H:%M')
            d.history = (getattr(d, 'history', '') + f'\n{stamp} V2 Q-80 add_site_crop.py: croptyp=[{cft}], '
                         f'pctcrop=[1.0] ({what})').strip()
        print('  written')
    if apply:
        os.makedirs(LOG, exist_ok=True)
        with open(f'{LOG}/md5_before_sitedata.txt', 'a') as f:
            f.write('\n'.join(md5) + '\n')
        print(f'md5 before in {LOG}/md5_before_sitedata.txt')
    else:
        print('dry run: nothing written (add --apply)')


if __name__ == '__main__':
    main()
