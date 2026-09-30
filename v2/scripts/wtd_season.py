"""Monthly climatology of the water level at the wetland towers: observed
(ledger D-2 table, m, + below surface, - ponded) against K4 (model level =
minus the pond depth above 1 mm, else the water-table depth). Q-18.
Usage: wtd_season.py [<case> <tag> [site ...]]; sites given replace the fen and bog groups."""
import sys
sys.path.insert(0, '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/scripts')
import numpy as np
from summarize_tag import OB, water_level
case, tag = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ('paper_v2/v260925r/sp_v2_c22', 'v2_k4_sasu5')
groups = {'fen': ['DE-Zrk', 'FI-Lom', 'FI-Sii', 'FR-LGt', 'SE-Deg', 'US-Los'],
          'bog': ['CA-SCB', 'DE-SfN', 'FI-Si2', 'JP-BBY', 'NZ-Kop']}
if len(sys.argv) > 3:
    groups = {'given': sys.argv[3:]}
for g, sites in groups.items():
    print(f'== {g}: month  obs  model  (median over towers of the monthly climatology, m)')
    O, M = [], []
    for s in sites:
        m = water_level(case, tag, s)
        j = m.merge(OB[OB.site == s][['year', 'month', 'wtd_obs']], on=['year', 'month'])
        o = j.groupby('month').wtd_obs.mean().reindex(range(1, 13))
        md = j.groupby('month').lev.mean().reindex(range(1, 13))
        O.append(o.values); M.append(md.values)
        print(f'   {s}: obs ' + ' '.join(f'{x:5.2f}' for x in o.values) + '\n          mod ' + ' '.join(f'{x:5.2f}' for x in md.values))
    O, M = np.array(O), np.array(M)
    print('   median obs ' + ' '.join(f'{x:5.2f}' for x in np.nanmedian(O, 0)))
    print('   median mod ' + ' '.join(f'{x:5.2f}' for x in np.nanmedian(M, 0)))
