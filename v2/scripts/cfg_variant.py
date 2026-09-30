#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Make a V2 site configuration variant in one call.
Copies v2/config/<base> to v2/config/<new> (unless <new> exists) and sets
KEY = VALUE in template.nml ('template') or ch4_parameter.nml ('ch4'): an
existing assignment of KEY is replaced in place, otherwise the comment and the
assignment are inserted before the '/' that closes the namelist group already
holding keys with KEY's prefix (DEF_METHANE% and the like; the first group for
a plain DEF_ key). Prints the diff against <base>.
Usage: cfg_variant.py <base> <new> template|ch4 <KEY> <VALUE> [comment]"""
import difflib
import os
import re
import shutil
import sys

CFG = '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/config'


def main():
    base, new, which, key, value = sys.argv[1:6]
    comment = sys.argv[6] if len(sys.argv) > 6 else ''
    src, dst = f'{CFG}/{base}', f'{CFG}/{new}'
    if not os.path.isdir(dst):
        shutil.copytree(src, dst)
    fname = {'template': 'template.nml', 'ch4': 'ch4_parameter.nml'}[which]
    path = f'{dst}/{fname}'
    lines = open(path).read().split('\n')
    pat = re.compile(r'^(\s*)' + re.escape(key) + r'\s*=')
    hit = [i for i, l in enumerate(lines) if pat.match(l)]
    indent = '    ' if which == 'ch4' else '   '
    if hit:
        i = hit[0]
        lines[i] = f'{pat.match(lines[i]).group(1)}{key} = {value}'
    else:
        prefix = key.split('%')[0] + '%' if '%' in key else ''
        start = next((i for i, l in enumerate(lines) if prefix and l.strip().startswith(prefix)), 0)
        end = next(i for i, l in enumerate(lines) if i >= start and l.strip() == '/')
        add = ([f'{indent}! {comment}'] if comment else []) + [f'{indent}{key} = {value}']
        lines[end:end] = add
    open(path, 'w').write('\n'.join(lines))
    old = open(f'{src}/{fname}').read().split('\n')
    sys.stdout.writelines(l + '\n' for l in difflib.unified_diff(old, lines, f'{base}/{fname}', f'{new}/{fname}', n=0, lineterm=''))


if __name__ == '__main__':
    main()
