"""All towers as one set (user 2026-09-26: the paper makes no main/validation
split, every change is a mechanism of the whole model): merges the per-site
lines of a main-set summary and a validation summary, and prints the medians
of KGE_ln, beta, alpha and r over all non-saline towers and by wetland type
(main-set classes and validation labels merged; salt marshes reported apart,
the model has no sulphate suppression).
Usage: all_sites.py <main summary> <validation summary> [label]"""
import re
import sys
from collections import defaultdict
import numpy as np

TYPE = {'bog': 'bog', 'Bog': 'bog', 'fen': 'fen', 'Fen': 'fen', 'marsh': 'marsh', 'Marsh': 'marsh',
        'trop_swamp': 'swamp', 'Swamp': 'swamp', 'Wet tundra': 'wet tundra', 'permafrost': 'permafrost',
        'Drained': 'drained', 'salt_marsh': 'salt marsh'}
PAT = re.compile(r'^([A-Z]{2}-[A-Za-z0-9]{3})\s+(.+?)\s+(-?[\d.]+|nan)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+|nan)\s*$')


def sites(path):
    out = {}
    for line in open(path):
        m = PAT.match(line.rstrip())
        if m and m.group(2).strip() in TYPE:
            out[m.group(1)] = (TYPE[m.group(2).strip()], *map(float, m.group(3, 4, 5, 6)))
    return out


def med(rows):
    a = np.array([r[1:] for r in rows], float)
    return ' / '.join(f'{x:5.2f}' for x in np.nanmedian(a, 0))


def main():
    s = {**sites(sys.argv[1]), **sites(sys.argv[2])}
    label = sys.argv[3] if len(sys.argv) > 3 else ''
    rows = [v for v in s.values() if v[0] != 'salt marsh' and np.isfinite(v[1])]
    print(f'# {label} all towers: KGE_ln / beta / alpha / r medians')
    print(f'all non-saline n={len(rows):2d}  {med(rows)}')
    g = defaultdict(list)
    for v in s.values():
        if np.isfinite(v[1]):
            g[v[0]].append(v)
    for t in ('bog', 'fen', 'marsh', 'wet tundra', 'swamp', 'drained', 'permafrost', 'salt marsh'):
        if g[t]:
            print(f'  {t:11s} n={len(g[t]):2d}  {med(g[t])}')


if __name__ == '__main__':
    main()
