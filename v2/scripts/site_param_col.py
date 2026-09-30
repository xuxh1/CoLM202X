#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Set one per-tower column of a V2 site configuration in one call.
Adds COLUMN to v2/config/<cfg>/<table> if missing (other towers left empty,
i.e. the namelist default) and writes VALUE for the listed towers; towers
not yet in the table get a row. Prints the towers set.
  table SITE_PARAMS.csv holds DEF_METHANE% keys (per-tower CH4 namelist),
  SITE_MAIN.csv holds keys of the main &nl_colm group (prep_tag.sh).
Usage: site_param_col.py <cfg> SITE_PARAMS.csv|SITE_MAIN.csv <COLUMN> <VALUE> <tower,tower,...>"""
import csv
import sys
from pathlib import Path

CFG = Path('/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/config')


def main():
    cfg, table, col, val, towers = sys.argv[1:6]
    p = CFG / cfg / table
    rows = list(csv.DictReader(open(p, newline=''))) if p.is_file() else []
    cols = list(rows[0].keys()) if rows else ['ID']
    if col not in cols:
        cols.append(col)
    want = [t.strip() for t in towers.split(',') if t.strip()]
    have = {r['ID'] for r in rows}
    rows += [{'ID': t} for t in want if t not in have]
    for r in rows:
        if r['ID'] in want:
            r[col] = val
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, cols, restval='')
        w.writeheader()
        w.writerows(rows)
    print(f'{p}: {col} = {val} for {", ".join(want)}')


if __name__ == '__main__':
    main()
