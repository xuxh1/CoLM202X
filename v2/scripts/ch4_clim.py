"""Monthly climatology of the CH4 flux at chosen towers for several site tags
against the tower (mg CH4 m-2 d-1), over the months the tower has, with the
score_series_v2 monthly pairing (V-43: does microtopography bring back the
spring peak of Q-4?).
Usage: ch4_clim.py <case:tag>[,<case:tag>...] <site> [site ...]"""
import sys
from collections import defaultdict
import numpy as np

sys.path.insert(0, '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/scripts')
import score_series_v2 as S                                          # noqa: E402

runs = [x.split(':') for x in sys.argv[1].split(',')]
for sid in sys.argv[2:]:
    obs = S.obs_monthly(sid)
    print(f'{sid}: month ' + ' '.join(f'{m:5d}' for m in range(1, 13)))
    ks = [k for k in obs if np.isfinite(obs[k])]
    rows = [('obs', {k: obs[k] for k in ks})]
    for case, tag in runs:
        mod = {k: v for k, v in
               S.S.mod_monthly(f'/share/home/dq076/mode/Methane/cases/{case}', tag, sid).items()}
        rows.append((tag, {k: mod[k] for k in ks if k in mod}))
    for name, d in rows:
        acc = defaultdict(list)
        for (y, m), v in d.items():
            acc[m].append(v)
        print(f'  {name[:16]:16s} ' + ' '.join(f'{np.mean(acc[m]):5.0f}' if acc[m] else '    .' for m in range(1, 13)))
