#!/share/home/dq076/software/miniconda3/envs/py311/bin/python
"""Site-lane scheduler of the V2 loop (v2/目标与判据.md section 5): keeps the
direct-run nodes busy from a queue of site tags instead of waiting for the
controlling session.

Queue v2/queue.tsv (tab separated, header line), one tag per row:
  id status tree tag cfg spin sites lib node started finished note
  status   pending | running | done | failed
  tree     case under cases/ (a built single-point tree)
  cfg      v2/config subdirectory (template.nml, ch4_parameter.nml,
           LIST_sites.csv, SITE_PARAMS.csv, SITE_MAIN.csv)
  spin     spin-up cycles written into the namelists
  sites    all, or a subset name from v2/config/SUBSETS.tsv
  lib      - or a frozen equilibrium (scripts/sites/freeze_spinup.py) to
           start stage 1 from (--spin-from; same commit as the tree)
Every POLL seconds: a running row whose sites_<tag>/_conf/STATUS.txt no
longer reads RUNNING is closed, scored (score77.py) and summarised into
v2/results/queue_summary.tsv; then pending rows are started on nodes with a
free slot (SLOTS tags per node, 1-min load below LOAD_MAX). Rows edited by
hand are picked up on the next poll. Stops when v2/queue.stop exists.
Usage: queue_runner.py   (run in the background on mgt01)"""
import csv
import datetime
import subprocess
import time
from pathlib import Path

V2 = Path('/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2')
CASES = Path('/share/home/dq076/mode/Methane/cases')
PROJ = Path('/share/home/dq076/mode/Methane')
Q = V2 / 'queue.tsv'
SUMMARY = V2 / 'results' / 'queue_summary.tsv'
STOP = V2 / 'queue.stop'
NODES = ('node110', 'node111', 'node112', 'node114')
SLOTS, LOAD_MAX, WORKERS, POLL = 2, 40, 24, 120
COLS = ('id', 'status', 'tree', 'tag', 'cfg', 'spin', 'sites', 'lib', 'node', 'started', 'finished', 'note')
CLASSES = ('all', 'wetland', 'bog', 'fen', 'marsh', 'swamp', 'wet tundra', 'drained', 'permafrost', 'rice')


def now():
    return datetime.datetime.now().strftime('%m-%d %H:%M')


def read_queue():
    with open(Q, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_queue(rows):
    # rows appended by hand while this pass held the queue in memory (scoring
    # takes minutes) would be lost; take them over before writing
    known = {r['id'] for r in rows}
    rows.extend(r for r in read_queue() if r['id'] not in known)
    tmp = Q.with_suffix('.tmp')
    with open(tmp, 'w', newline='') as f:
        w = csv.DictWriter(f, COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    tmp.replace(Q)


def subsets():
    out = {}
    p = V2 / 'config' / 'SUBSETS.tsv'
    for line in p.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            name, sites = line.split('\t', 1)
            out[name.strip()] = sites.strip()
    return out


def status_of(r):
    s = CASES / r['tree'] / f"sites_{r['tag']}" / '_conf' / 'STATUS.txt'
    return s.read_text().split('\n', 1)[0].strip() if s.is_file() else ''


def load(node):
    try:
        out = subprocess.run(['ssh', '-o', 'ConnectTimeout=5', node, 'cat', '/proc/loadavg'],
                             capture_output=True, text=True, timeout=20).stdout
        return float(out.split()[0])
    except Exception:
        return 99.


def score(r):
    p = V2 / 'results' / f"scores77_{r['tag']}.txt"
    # node114 alone failed silently when overloaded (v2_p9w, v2_p9ph had no
    # score file): fall back to the node that ran the tag, and log a failure
    for node in dict.fromkeys(['node114', r.get('node') or 'node114']):
        res = subprocess.run([str(V2 / 'scripts' / 'run_on_node.sh'), node, 'score77.py', r['tree'], r['tag']],
                             cwd=V2 / 'scripts', capture_output=True, text=True, timeout=3600)
        if p.is_file() and p.stat().st_mtime >= time.time() - 3600:
            break
        print(f"{datetime.datetime.now():%m-%d %H:%M} score {r['tag']} on {node} left no fresh {p.name}: {res.stderr.strip()[-300:]}", flush=True)
    med = {}
    if p.is_file():
        for line in p.read_text().splitlines():
            for c in CLASSES:
                if line.startswith(c + ' ') and 'median' not in line:
                    parts = line[len(c):].split()
                    if len(parts) >= 5:
                        med[c] = f'{parts[1]}/{parts[2]}/{parts[3]}/{parts[4]}'
    new = not SUMMARY.is_file()
    with open(SUMMARY, 'a') as f:
        if new:
            f.write('tag\ttree\tcfg\tsites\tfinished\t' + '\t'.join(CLASSES) + '\n')
        f.write('\t'.join([r['tag'], r['tree'], r['cfg'], r['sites'], r['finished']]
                          + [med.get(c, '') for c in CLASSES]) + '\n')


def start(r, node, subs):
    cfg = V2 / 'config' / r['cfg']
    prep = subprocess.run([str(V2 / 'scripts' / 'prep_tag.sh'), r['tree'], r['tag'], r['cfg'],
                           str(cfg / 'LIST_sites.csv'), r['spin'], str(cfg / 'SITE_PARAMS.csv'),
                           str(cfg / 'SITE_MAIN.csv')], capture_output=True, text=True, timeout=600)
    if prep.returncode != 0:
        return f'prep failed: {prep.stderr[-200:]}'
    extra = []
    if r['sites'] not in ('', 'all'):
        extra += ['--sites', subs[r['sites']]]
    if r['lib'] not in ('', '-'):
        extra += ['--spin-from', r['lib']]
    env = f'LOAD_MAX={LOAD_MAX} WORKERS={WORKERS}'
    cmd = f"cd {PROJ} && {env} scripts/sites/run_sites_direct.sh {node} {r['tree']} {r['tag']} {' '.join(extra)}"
    run = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True, timeout=600)
    if 'log:' not in run.stdout:
        return f'launch failed: {(run.stdout + run.stderr)[-200:]}'
    return ''


def main():
    while not STOP.exists():
        rows = read_queue()
        subs = subsets()
        changed = False
        for r in rows:
            if r['status'] == 'running':
                st = status_of(r)
                if st and st != 'RUNNING':
                    r['status'] = 'done' if st == 'DONE' else 'failed'
                    r['finished'] = now()
                    changed = True
                    write_queue(rows)
                    score(r)
        busy = {n: sum(1 for r in rows if r['status'] == 'running' and r['node'] == n) for n in NODES}
        for r in rows:
            if r['status'] != 'pending':
                continue
            free = [n for n in NODES if busy[n] < SLOTS]
            free = [n for n in free if load(n) < LOAD_MAX]
            if not free:
                break
            node = min(free, key=lambda n: busy[n])
            err = start(r, node, subs)
            if err:
                r['status'] = 'failed'
                r['note'] = (r['note'] + ' | ' + err).strip(' |')
            else:
                r['status'], r['node'], r['started'] = 'running', node, now()
                busy[node] += 1
            changed = True
            write_queue(rows)
        if changed:
            write_queue(rows)
        time.sleep(POLL)


if __name__ == '__main__':
    main()
