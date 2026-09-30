#!/bin/bash
# Write a site tag's namelists from a V2 config directory.
# Usage: prep_tag.sh <case under cases/> <tag> <v2/config subdir> <sitelist csv> <spinup> [site-param csv] [site-main csv]
#   site-param csv: per-site DEF_METHANE% overrides (gen_site_nmls --site-param-csv, CH4 namelist).
#   site-main csv: per-site keys for the main &nl_colm group (ID column plus one column per key,
#   e.g. DEF_WETLAND_PEAT_SHARE_SITE); gen_site_nmls only knows SITE_ph/salinity/wetland_class,
#   so these are written into both of the tower's namelists here, after generation.
set -eo pipefail
CASE=${1:?case}; TAG=${2:?tag}; CFG=${3:?config}; LIST=${4:?sitelist}; SPIN=${5:?spinup}; SP=${6:-}; SM=${7:-}
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
PY=/share/home/dq076/software/miniconda3/envs/py311/bin/python
cd /share/home/dq076/mode/Methane
EXTRA=(); [ -n "$SP" ] && EXTRA=(--site-param-csv "$SP")
$PY scripts/sites/gen_site_nmls.py --case "$CASE" --tag "$TAG" \
  --sitelist "$LIST" --template "$V2/config/$CFG/template.nml" --ch4param "$V2/config/$CFG/ch4_parameter.nml" \
  --spinup "$SPIN" "${EXTRA[@]}" 2>&1 | tail -2
if [ -n "$SM" ]; then
  D=cases/$CASE/sites_$TAG/_conf
  cp "$SM" "$D/SITE_MAIN.csv"
  $PY - "$SM" "$D/site_nml" <<'PYEOF'
import csv, re, sys, pathlib
sm, nd = sys.argv[1], pathlib.Path(sys.argv[2])
n = 0
for row in csv.DictReader(open(sm, newline='')):
    sid = row.pop('ID').strip()
    kv = [(k.strip(), v.strip()) for k, v in row.items() if v and v.strip()]
    for f in (nd / f'{sid}.nml', nd / f'{sid}_spin.nml'):
        if f.exists():
            t = f.read_text()
            assert t.startswith('&nl_colm\n'), f
            # a key the template already sets is replaced in place (a second
            # assignment above it would be overridden by the template's);
            # other keys go right after &nl_colm
            new = ''
            for k, v in kv:
                pat = re.compile(r'^([ \t]*)' + re.escape(k) + r'[ \t]*=.*$', re.M | re.I)
                if pat.search(t):
                    t = pat.sub(lambda m: f'{m.group(1)}{k} = {v}', t, count=1)
                else:
                    new += f'   {k} = {v}\n'
            f.write_text(t.replace('&nl_colm\n', '&nl_colm\n' + new, 1))
            n += 1
print(f'  site main keys written into {n} namelists from {sm}')
PYEOF
fi
# Crop towers: CN crop phenology sets gdd020 only on the first step of 1 January
# (CNPhenology.F90:318-337), so a spin-up window that stops at 1 January 00:00
# (a forcing that starts on 2 January of a leap year) never initialises it and
# rice matures in 50-68 days. Configs carrying SPIN_CROSS_JAN1 move such a
# window end to 2 January.
if [ -f "$V2/config/$CFG/SPIN_CROSS_JAN1" ]; then
  $PY - "cases/$CASE/sites_$TAG/_conf/site_nml" <<'PYEOF'
import re, sys, pathlib
n = 0
for f in pathlib.Path(sys.argv[1]).glob('*_spin.nml'):
    t = f.read_text()
    mo = re.search(r'DEF_simulation_time%spinup_month\s*=\s*(\d+)', t)
    dy = re.search(r'DEF_simulation_time%spinup_day\s*=\s*(\d+)', t)
    if mo and dy and int(mo.group(1)) == 1 and int(dy.group(1)) == 1:
        f.write_text(re.sub(r'(DEF_simulation_time%spinup_day\s*=\s*)1\b', r'\g<1>2', t)); n += 1
        print(f'  spin window of {f.stem} moved to 2 January')
print(f'  spin windows checked for 1 January: {n} moved')
PYEOF
fi
