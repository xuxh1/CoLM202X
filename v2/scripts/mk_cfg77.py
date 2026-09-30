#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Build a V2 site configuration for all 77 towers of scripts/sites/LIST_sites_77.csv
from a 46-tower configuration (its own LIST/SITE_PARAMS/SITE_MAIN for the former
main set plus the shared validation tables), for the model-optimisation loop
(user 2026-09-27: score all towers, not only the 23 with water table).

- LIST_sites.csv: the 77 rows; SITE_wetland_class, SITE_ph, SITE_salinity blank,
  as in every V2 list (V2 reads the class axes from SITE_PARAMS / SITE_MAIN).
- SITE_PARAMS.csv: former main rows, then validation rows, then the three salt
  marshes not in either set (US-EDN, US-LA1, US-MRM: no trees, marsh LAI cap 3.0,
  as the other salt marshes). Non-wetland towers take the template values.
- SITE_MAIN.csv: peat share and USE_SITE_LAI as before (salt marshes mineral);
  USE_SITE_pctcrop = .true. for the seven rice towers (single irrigated-rice
  patch from croptyp/pctcrop in the master site files, D-24).
- ch4_parameter.nml: the base file plus the paper rice keys C-24 and C-25b
  when absent (they act on rice columns only).
- SPIN_CROSS_JAN1 marker so rice spin windows cross 1 January (prep_tag.sh).
Usage: mk_cfg77.py <base config> <new config>"""
import os
import re
import shutil
import sys

import pandas as pd

CFG = '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/config'
L77 = '/share/home/dq076/mode/Methane/scripts/sites/LIST_sites_77.csv'
SALT_NEW = ('US-EDN', 'US-LA1', 'US-MRM')
RICE_KEYS = [('DEF_METHANE%rice_aere_override', '.false.',
              '! C-24 (paper V2): paddy rice keeps the CLM4Me aerenchyma geometry (Riley et al. 2011).'),
             ('DEF_METHANE%rice_aereoxid', '0.65',
              '! C-25b (paper V2): rhizosphere oxidation of the aerenchyma CH4 in paddy rice.')]


def main():
    base, new = sys.argv[1:3]
    src, dst = f'{CFG}/{base}', f'{CFG}/{new}'
    if os.path.exists(dst):
        sys.exit(f'{dst} exists')
    os.makedirs(dst)
    shutil.copy2(f'{src}/template.nml', dst)
    a = pd.read_csv(L77, dtype=str)
    for c in ('SITE_wetland_class', 'SITE_ph', 'SITE_salinity'):
        a[c] = ''
    a.to_csv(f'{dst}/LIST_sites.csv', index=False)

    pm = pd.read_csv(f'{src}/SITE_PARAMS.csv', dtype=str)
    pv_path = f'{src}/SITE_PARAMS_val23.csv'
    pv = pd.read_csv(pv_path if os.path.exists(pv_path) else f'{CFG}/SITE_PARAMS_val23_c13.csv', dtype=str)
    salt = pd.DataFrame({'ID': SALT_NEW, pm.columns[1]: '0.0', pm.columns[2]: '3.0'})
    p = pd.concat([pm, pv, salt], ignore_index=True)
    assert not p.ID.duplicated().any()
    p.to_csv(f'{dst}/SITE_PARAMS.csv', index=False)

    mm = pd.read_csv(f'{src}/SITE_MAIN.csv', dtype=str)
    mv = pd.read_csv(f'{CFG}/SITE_MAIN_val23.csv', dtype=str)
    m = pd.concat([mm, mv, pd.DataFrame({'ID': SALT_NEW, 'DEF_WETLAND_PEAT_SHARE_SITE': '0.0'})],
                  ignore_index=True)
    rice = a.loc[a.SITE_CLASSIFICATION == 'Rice', 'ID'].tolist()
    m = pd.concat([m, pd.DataFrame({'ID': rice})], ignore_index=True)
    m['USE_SITE_pctcrop'] = ''
    m.loc[m.ID.isin(rice), 'USE_SITE_pctcrop'] = '.true.'
    assert not m.ID.duplicated().any()
    m.to_csv(f'{dst}/SITE_MAIN.csv', index=False)

    t = open(f'{src}/ch4_parameter.nml').read()
    add = ''
    for k, v, c in RICE_KEYS:
        if not re.search(r'^\s*' + re.escape(k) + r'\s*=', t, re.M):
            add += f'    {c}\n    {k} = {v}\n'
    if add:
        # before the '/' that closes &nl_colm_methane_parameter (not the file's last group)
        g = t.index('&nl_colm_methane_parameter')
        i = g + re.search(r'^\s*/\s*$', t[g:], re.M).start()
        t = t[:i] + add + t[i:]
    open(f'{dst}/ch4_parameter.nml', 'w').write(t)
    open(f'{dst}/SPIN_CROSS_JAN1', 'w').write('rice spin windows cross 1 January (prep_tag.sh)\n')
    print(f'{dst}: {len(a)} towers, SITE_PARAMS {len(p)}, SITE_MAIN {len(m)} (rice {len(rice)})')


if __name__ == '__main__':
    main()
