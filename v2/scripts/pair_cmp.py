#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Compare a site tag with a baseline tag on the towers both scored, so a
subset screen is judged against the baseline on the same towers (the class
medians in queue_summary.tsv are over each tag's own tower set).
Reads v2/results/scores77_<tag>.txt of both tags.
Prints per tower KGE_ln base -> new with beta and r, then per class the
medians of the four scores for both tags on the common towers.
Usage: pair_cmp.py <base tag> <new tag> [min |dKGE| to list, default 0.05]"""
import sys
from pathlib import Path

import numpy as np

R = Path('/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/results')


def read(tag):
    out = {}
    for line in (R / f'scores77_{tag}.txt').read_text().splitlines():
        p = line.split()
        if len(p) < 9 or line.startswith(('#', 'site', 'group')):
            continue
        try:
            vals = [float(x) for x in p[-7:]]
        except ValueError:
            continue
        # class names may hold a space (wet tundra)
        out[p[0]] = (' '.join(p[1:-7]), vals[1:5])
    return out


def main():
    base, new = read(sys.argv[1]), read(sys.argv[2])
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    common = sorted(set(base) & set(new))
    print(f'{len(common)} common towers; listed where |dKGE_ln| >= {thr}')
    print(f"{'site':8s} {'class':12s} {'KGE base':>8s} {'new':>6s} {'beta b':>7s} {'new':>6s} {'r b':>6s} {'new':>6s}")
    for s in sorted(common, key=lambda s: new[s][1][0] - base[s][1][0]):
        b, n = base[s][1], new[s][1]
        if abs(n[0] - b[0]) >= thr:
            print(f'{s:8s} {base[s][0]:12s} {b[0]:8.2f} {n[0]:6.2f} {b[1]:7.2f} {n[1]:6.2f} {b[3]:6.2f} {n[3]:6.2f}')
    classes = sorted({base[s][0] for s in common})
    print(f"\n{'class':12s} {'n':>3s}  base KGE/beta/alpha/r        new KGE/beta/alpha/r   up/down")
    for c in ['all'] + classes:
        ss = [s for s in common if c == 'all' or base[s][0] == c]
        mb = np.median([base[s][1] for s in ss], axis=0)
        mn = np.median([new[s][1] for s in ss], axis=0)
        up = sum(new[s][1][0] > base[s][1][0] + 0.02 for s in ss)
        dn = sum(new[s][1][0] < base[s][1][0] - 0.02 for s in ss)
        print(f"{c:12s} {len(ss):3d}  {'/'.join(f'{x:.2f}' for x in mb):26s} {'/'.join(f'{x:.2f}' for x in mn):26s} {up}/{dn}")


if __name__ == '__main__':
    main()
