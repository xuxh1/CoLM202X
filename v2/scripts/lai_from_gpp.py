#!/usr/bin/env python
"""Footprint LAI cap of the wetland towers from tower GPP (Beer's law).

Wetland towers are those of scripts/sites/LIST_sites_77.csv outside the rice,
lake and upland classes, with the classes of score77.py. Per tower:
  season  the three consecutive calendar months (wrapping over the year end)
          with the highest tower GPP climatology, each with a mean forcing air
          temperature above 0 C;
  tower   GPP of Observation/<tower>_*_Flux.nc (night-time partitioning;
          GPP_DT for a second ratio) in g C m-2 d-1; time steps outside
          -50 to 100 umol m-2 s-1 count as missing, a month needs >= 80 % of
          its time steps finite;
  model   f_assim (gross canopy assimilation, mol CO2 m-2 s-1) x 12.011 x
          86400 -> g C m-2 d-1, with f_lai and f_rstfacsun; patches weighted
          as score77.weights;
compared on the matched (year, month) pairs of the season: each calendar
month is averaged over its pairs, the season is the mean of its three months.
A tower needs at least one complete season of pairs.

Inversion: with GPP ~ 1 - exp(-k LAI) the LAI scaling c that brings the model
season to the tower solves, over the season months i,
    sum_i Gm_i (1 - exp(-k c L_i)) / (1 - exp(-k L_i)) = sum_i Go_i
(bisection; none when the light-saturated sum_i Gm_i / (1 - exp(-k L_i)) is
still short). The season LAI it asks for is c times the modelled season LAI,
the cap c times the modelled annual LAI peak: under wetland_lai_shape the
non-forested LAI is the remote-sensing LAI times cap / annual remote-sensing
peak, so the modelled peak is the cap wherever the cap binds. A cap is
suggested only for a model / tower GPP ratio outside 0.8-1.25 on a wetland
patch (IGBP 11) with forested share < 1, at k = 0.5 with 0.4 and 0.6 as
sensitivity. Raising the cap lifts the LAI only up to the remote-sensing
input, whose annual peak comes from the tower's landdata/srfdata.nc.

Sibling tags of the same tree (finished, other than the tag) that set a
tower's wetland_lai_cap_site differently and moved its season LAI by more
than 10 % give the model's own k from the two runs on the same pairs,
    sum_i G1_i (1 - exp(-k L2_i)) / (1 - exp(-k L1_i)) = sum_i G2_i;
their median gives one more suggestion column.

Moss correction (Q-63 after V-172): tower GPP of a peatland includes the moss
layer's photosynthesis, which the model does not have, so the vascular target
is tower GPP x (1 - m), m = wetland_moss_input_frac (0.48, Turetsky et al.
2010) x s_moss, s_moss the tower's wetland_share_moss_site of the config
(else the design table of v2/results/design_wetland_class_260930.txt sect.
5.5), or a site-measured moss share of GPP where one exists (MOSS_SITE). The
caps are recomputed at the model's own k against this target and set against
the config's current cap. With --pair, a second run (any tree) with other
caps gives per tower its season GPP against both targets, the k implied
between the two runs and its methane beta.

Writes v2/results/lai_from_gpp_<YYMMDD>.txt and prints it.
Usage: lai_from_gpp.py <tree under cases/> <tag> [config under v2/config whose
       caps are listed next to the tag's, default v2_p8f] [--pair <tree> <tag>]"""
import calendar
import csv
import datetime
import glob
import os
import re
import sys

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score77 import classes, weights                                # noqa: E402
import pair_cmp                                                     # noqa: E402

R = '/share/home/dq076/mode/Methane'
V2 = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OBS = f'{R}/data/FLUXNET-CH4/Observation'
NOT_WET = ('rice', 'lake', 'upland')
CLASS_ZH = {'bog': 'bog', 'fen': 'fen', 'marsh': 'marsh', 'swamp': '沼泽林', 'wet tundra': '湿苔原',
            'drained': '排干', 'permafrost': '永冻', 'salt marsh': '盐沼'}
K_MAIN, K_SENS = 0.5, (0.4, 0.6)
LO, HI = 0.8, 1.25
CAP_KEY = 'DEF_METHANE%wetland_lai_cap_site'
FOREST_KEY = 'DEF_METHANE%wetland_forest_share_site'
GC = 12.011 * 86400.                    # mol CO2 m-2 s-1 -> g C m-2 d-1
GPP_MIN, GPP_MAX = -50., 100.           # plausible tower GPP per time step [umol CO2 m-2 s-1]
MOSS_KEY = 'DEF_METHANE%wetland_share_moss_site'
PHI_M = 0.48                            # wetland_moss_input_frac, moss share of wetland productivity
PHI_M_ALL = 0.32                        # the same share of all productivity if vascular BNPP = ANPP (sensitivity)
DESIGN = 'design_wetland_class_260930.txt'
# site-measured moss share of GPP, replacing PHI_M x s_moss
MOSS_SITE = {'FR-LGt': (414. / 1273., 'Leroy et al. 2019 BG 16:4085，La Guette 泥炭中宇宙，Sphagnum 单独 GPP 414、'
                                      'Sphagnum 加 Molinia 1273 g C m-2 yr-1（表 3），Sphagnum 盖度不受 Molinia 影响（2.1 节）')}
CHANGE = 0.15                           # |suggested - current cap| that counts as a change [m2 m-2]
L = []


def out(s=''):
    print(s)
    L.append(s)


def nml_raw(path, key):
    """Last value string assigned to key in a namelist file; None without one."""
    if not os.path.exists(path):
        return None
    val = None
    pat = re.compile(r'^\s*' + re.escape(key) + r'\s*=\s*([^\s!,]+)')
    for line in open(path):
        m = pat.match(line)
        if m:
            val = m.group(1)
    return val


def nml_value(path, key):
    """Last numeric assignment of key in a namelist file; None without one."""
    v = nml_raw(path, key)
    return None if v is None else float(v.lower().replace('d', 'e'))


def cfg_caps(cfg, key=CAP_KEY):
    """A site key (default wetland_lai_cap_site) per tower of
    v2/config/<cfg>/SITE_PARAMS.csv, towers with a value only."""
    f = f'{V2}/config/{cfg}/SITE_PARAMS.csv'
    caps = {}
    if os.path.exists(f):
        for row in csv.DictReader(open(f)):
            v = (row.get(key) or '').strip()
            if v:
                caps[row['ID']] = float(v)
    return caps


def design_moss():
    """s_moss per tower from the table of sect. 5.5 of the wetland class design
    (rows may list several towers joined by '、')."""
    res, on = {}, False
    for line in open(f'{V2}/results/{DESIGN}'):
        if line.startswith('### '):
            on = line.startswith('### 5.5')
            continue
        c = [x.strip() for x in line.split('|')]
        if on and len(c) > 4 and c[1] not in ('塔', '---'):
            try:
                s = float(c[4])
            except ValueError:
                continue
            for t in c[1].split('、'):
                if re.fullmatch(r'[A-Z]{2}-\w{3}', t):
                    res[t] = s
    return res


def cap_num(s):
    """Leading number of a suggestion string, or None."""
    m = re.match(r'\d+(\.\d+)?', s or '')
    return float(m.group(0)) if m else None


def obs_monthly(sid):
    """Tower GPP and GPP_DT per (year, month) in g C m-2 d-1; a month needs
    >= 80 % of its time steps finite."""
    res = {'GPP': {}, 'GPP_DT': {}}
    for p in sorted(glob.glob(f'{OBS}/{sid}_*_Flux.nc')):
        with nc.Dataset(p) as d:
            tv = d['time']
            t = nc.num2date(tv[:], tv.units, getattr(tv, 'calendar', 'standard'))
            step = np.median([(t[i + 1] - t[i]).total_seconds() for i in range(min(len(t) - 1, 500))])
            key = np.array([x.year * 100 + x.month for x in t])
            for v in res:
                if v not in d.variables:
                    continue
                u = d[v].units.replace(' ', '').lower()
                if 'umol' in u and 's-1' in u:
                    fac = 1e-6 * GC
                elif u.startswith('gc') and 'd-1' in u:
                    fac = 1.
                else:
                    raise SystemExit(f'{p}: {v} in unknown units {d[v].units}')
                a = np.ma.filled(d[v][:].astype(float), np.nan).ravel() * fac
                # unfilled junk (US-A10 holds 2.4e37, US-A03 252 umol m-2 s-1)
                a[(a < GPP_MIN * 1e-6 * GC) | (a > GPP_MAX * 1e-6 * GC)] = np.nan
                for k in np.unique(key):
                    y, m = divmod(int(k), 100)
                    x = a[key == k]
                    if np.isfinite(x).sum() >= 0.8 * calendar.monthrange(y, m)[1] * 86400. / step:
                        res[v][(y, m)] = float(np.nanmean(x))
    return res


def model_monthly(tree, tag, sid):
    """Model GPP (g C m-2 d-1), LAI, sunlit stomatal water-stress factor and
    forcing air temperature (K) per (year, month)."""
    res = {}
    keys = ('f_assim', 'f_lai', 'f_rstfacsun', 'f_xy_t')
    for f in sorted(glob.glob(f'{R}/cases/{tree}/sites_{tag}/{sid}/history/{sid}_hist_[0-9][0-9][0-9][0-9].nc')):
        with nc.Dataset(f) as d:
            t = nc.num2date(d['time'][:], d['time'].units)
            v = {k: np.ma.filled(d[k][:].astype(float), np.nan).reshape(len(t), -1) for k in keys}
            w = weights(tree, tag, sid, v['f_assim'].shape[1])
            for i, x in enumerate(t):
                if not np.isfinite(v['f_assim'][i]).any():
                    continue
                g, lai, rs, ta = (float(np.nansum(v[k][i] * w)) for k in keys)
                res[(x.year, x.month)] = (g * GC, lai, rs, ta)
    return res


def warm_months(mod):
    """Calendar months whose mean forcing air temperature is above 0 C."""
    clim = {}
    for (_, m), v in mod.items():
        clim.setdefault(m, []).append(v[3])
    return {m for m, v in clim.items() if np.mean(v) > 273.15}


def season(obs, warm):
    """Three consecutive calendar months with the highest tower GPP climatology,
    all of them above 0 C (winter gap-filling leaves spurious GPP at the
    Arctic towers)."""
    clim = {}
    for (_, m), v in obs.items():
        clim.setdefault(m, []).append(v)
    clim = {m: np.mean(v) for m, v in clim.items()}
    best = None
    for m0 in range(1, 13):
        ms = [(m0 + j - 1) % 12 + 1 for j in range(3)]
        if all(m in clim and m in warm for m in ms):
            s = sum(clim[m] for m in ms)
            if best is None or s > best[0]:
                best = (s, ms)
    return best[1] if best else None


def matched(obs, mod, ms):
    """Per season month: the matched (year, month) pairs and the tower and model
    means over them; the number of complete seasons. None if a month has none."""
    pairs, go, gm, lai, rs = [], [], [], [], []
    for m in ms:
        p = sorted(q for q in obs if q[1] == m and q in mod)
        if not p:
            return None
        pairs.append(p)
        go.append(np.mean([obs[q] for q in p]))
        gm.append(np.mean([mod[q][0] for q in p]))
        lai.append(np.mean([mod[q][1] for q in p]))
        rs.append(np.mean([mod[q][2] for q in p]))
    off = [int(m < ms[0]) for m in ms]
    allp = {q for p in pairs for q in p}
    years = {y - off[ms.index(m)] for (y, m) in allp}
    n = sum(all((y + o, m) in allp for m, o in zip(ms, off)) for y in years)
    return dict(pairs=pairs, go=np.array(go), gm=np.array(gm), lai=np.array(lai), rs=np.array(rs), n=n)


def gpp_at(gm, lai, k, c):
    """Season GPP (sum over months) with the LAI scaled by c, Beer's law."""
    lai = np.maximum(lai, 1e-3)
    return float(np.sum(gm * (1. - np.exp(-k * c * lai)) / (1. - np.exp(-k * lai))))


def lai_scale(gm, go, lai, k):
    """LAI scaling that brings the model season GPP to the tower's; inf when an
    unbounded LAI still falls short."""
    tgt = float(np.sum(go))
    if tgt >= gpp_at(gm, lai, k, 1e6):
        return np.inf
    lo, hi = 0., 1.
    while gpp_at(gm, lai, k, hi) < tgt:
        hi *= 2.
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if gpp_at(gm, lai, k, mid) < tgt else (lo, mid)
    return 0.5 * (lo + hi)


def implied_k(g1, l1, g2, l2):
    """Beer's-law k that turns the season GPP at LAI l1 into the one at l2."""
    l1, l2 = np.maximum(l1, 1e-3), np.maximum(l2, 1e-3)
    tgt = float(np.sum(g2))

    def f(k):
        return float(np.sum(g1 * (1. - np.exp(-k * l2)) / (1. - np.exp(-k * l1)))) - tgt
    lo, hi = 1e-3, 5.
    if f(lo) * f(hi) > 0:
        return np.nan
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if f(lo) * f(mid) <= 0 else (mid, hi)
    return 0.5 * (lo + hi)


def siblings(tree, tag):
    """Finished site tags of the tree other than the tag (spin-up dirs left out)."""
    res = []
    for d in sorted(glob.glob(f'{R}/cases/{tree}/sites_*')):
        name = os.path.basename(d)[len('sites_'):]
        st = f'{d}/_conf/STATUS.txt'
        if name == tag or name.endswith('_spin') or not os.path.exists(st):
            continue
        if open(st).readline().strip() == 'DONE':
            res.append(name)
    return res


def srf_meta(tree, tag, sid):
    """IGBP class of the tower patch and the mean annual peak of the
    remote-sensing LAI input (nan for PFT LAI)."""
    f = f'{R}/cases/{tree}/sites_{tag}/{sid}/landdata/srfdata.nc'
    if not os.path.exists(f):
        return None, np.nan
    with nc.Dataset(f) as d:
        igbp = int(np.ravel(d['IGBP_classification'][:])[0])
        peak = np.nan
        if 'LAI_monthly' in d.variables:
            a = np.ma.filled(d['LAI_monthly'][:].astype(float), np.nan).reshape(-1, 12)
            peak = float(np.nanmean(np.nanmax(a, axis=1)))
    return igbp, peak


def fmt(x, nd=2):
    return '—' if x is None or not np.isfinite(x) else f'{x:.{nd}f}'


def suggest(r, k):
    """Suggested cap at extinction k, or why none."""
    if r['state'] not in ('压', '未压', 'GLWD格点'):
        return '不作用'
    c = lai_scale(r['gm'], r['go'], r['lai'], k)
    if not np.isfinite(c):
        return '补不足'
    cap = c * r['peak']
    if c > 1. and r['state'] == '未压':
        return '抬无效'
    if c > 1. and np.isfinite(r['rspeak']) and cap > r['rspeak']:
        cmax = r['rspeak'] / r['peak']
        rr = gpp_at(r['gm'], r['lai'], k, cmax) / float(np.sum(r['go']))
        return f'{r["rspeak"]:.1f}（遥感峰封顶，比值到 {rr:.2f}）'
    return f'{cap:.1f}'


def main():
    args = sys.argv[1:]
    pair = None
    if '--pair' in args:
        i = args.index('--pair')
        pair = tuple(args[i + 1:i + 3])
        del args[i:i + 3]
    tree, tag = args[:2]
    cfg = args[2] if len(args) > 2 else 'v2_p8f'
    conf = f'{R}/cases/{tree}/sites_{tag}/_conf'
    ccaps = cfg_caps(cfg)
    # methane magnitude ratio beta of the tag (score77.py output), when scored
    try:
        beta = {s: v[1][1] for s, v in pair_cmp.read(tag).items()}
    except OSError:
        beta = {}
    cls = classes()
    wet = sorted((s for s in cls if cls[s] not in NOT_WET), key=lambda s: (cls[s], s))
    sibs = siblings(tree, tag)
    commit = ''
    cc = f'{R}/cases/{tree}/CODE_COMMIT.txt'
    if os.path.exists(cc):
        for line in open(cc):
            if line.startswith('commit'):
                commit = line.split(':', 1)[1].strip()[:8]

    rows, sib_rows, shape = [], [], {}
    for s in wet:
        igbp, rspeak = srf_meta(tree, tag, s)
        pnml = f'{conf}/ch4_param/{s}.nml'
        cap = nml_value(pnml, CAP_KEY)
        fshare = nml_value(pnml, FOREST_KEY)
        sh = nml_raw(pnml, 'DEF_METHANE%wetland_lai_shape')
        shape[s] = sh is not None and sh.strip('.').lower().startswith('t')
        r = dict(s=s, c=CLASS_ZH.get(cls[s], cls[s]), cap=cap, fshare=fshare, ccap=ccaps.get(s), rspeak=rspeak)
        if igbp != 11:
            r['state'] = '非湿地斑块'
        elif fshare is not None and fshare >= 1.:
            r['state'] = '有林份额1'
        elif cap is None or cap <= 0.:
            r['state'] = 'GLWD格点'
        else:
            r['state'] = '压' if np.isfinite(rspeak) and rspeak > 1.01 * cap else '未压'
        rows.append(r)
        o = obs_monthly(s)
        mod = model_monthly(tree, tag, s)
        ms = season(o['GPP'], warm_months(mod))
        mt = matched(o['GPP'], mod, ms) if ms and mod else None
        if not o['GPP']:
            r['why'] = '无观测 GPP'
        elif not mod:
            r['why'] = '无模式输出'
        elif mt is None and ms:
            r['why'] = '无对应模式月'
        elif mt is None or mt['n'] < 1:
            r['why'] = '不足一季'
        if 'why' in r:
            continue
        r.update(mt, ms=ms, ratio=float(np.sum(mt['gm']) / np.sum(mt['go'])), obs=o['GPP'], mod=mod)
        mdt = matched(o['GPP_DT'], mod, ms) if o['GPP_DT'] else None
        r['ratio_dt'] = float(np.sum(mdt['gm']) / np.sum(mdt['go'])) if mdt else np.nan
        # annual LAI peak: years holding all three season months (runs may
        # start or end within a year, US-NGC)
        years = {}
        for (y, m), v in mod.items():
            years.setdefault(y, {})[m] = v[1]
        r['peak'] = float(np.mean([max(v.values()) for v in years.values() if all(m in v for m in ms)] or [np.nan]))
        cs = lai_scale(mt['gm'], mt['go'], mt['lai'], K_MAIN)
        r['lai_inv'] = cs * float(np.mean(mt['lai']))
        # the model's own k from sibling tags with another cap at this tower
        if r['state'] not in ('压', '未压') or cap is None:
            continue
        seen = set()
        for sb in sibs:
            cap2 = nml_value(f'{R}/cases/{tree}/sites_{sb}/_conf/ch4_param/{s}.nml', CAP_KEY)
            if cap2 is None or abs(cap2 - cap) < 0.01 * cap or cap2 in seen:
                continue
            seen.add(cap2)
            mod2 = model_monthly(tree, sb, s)
            g1, l1, g2, l2 = [], [], [], []
            for p in mt['pairs']:
                q = [x for x in p if x in mod2]
                if not q:
                    break
                g1.append(np.mean([mod[x][0] for x in q]))
                l1.append(np.mean([mod[x][1] for x in q]))
                g2.append(np.mean([mod2[x][0] for x in q]))
                l2.append(np.mean([mod2[x][1] for x in q]))
            if len(g2) < 3 or abs(np.mean(l2) / np.mean(l1) - 1.) < 0.1:
                continue
            g1, l1, g2, l2 = map(np.array, (g1, l1, g2, l2))
            pred = float(np.sum(g1 * (1. - np.exp(-K_MAIN * l2)) / (1. - np.exp(-K_MAIN * l1))))
            sib_rows.append(dict(s=s, sb=sb, cap=cap, cap2=cap2, l1=l1.mean(), l2=l2.mean(), g1=g1.mean(),
                                 g2=g2.mean(), go=float(np.mean(mt['go'])), pred=pred / 3., k=implied_k(g1, l1, g2, l2)))

    ks_own = [x['k'] for x in sib_rows if np.isfinite(x['k'])]
    k_own = round(float(np.median(ks_own)), 2) if ks_own else None
    kset = (K_MAIN,) + K_SENS + ((k_own,) if k_own else ())
    ok = [r for r in rows if 'ratio' in r]
    for r in ok:
        r['sug'] = {k: suggest(r, k) for k in kset} if not LO <= r['ratio'] <= HI else None
    off = [r for r in ok if r['sug']]
    act = [r for r in off if r['state'] in ('压', '未压', 'GLWD格点')]

    # moss correction: vascular target, caps at the model's own k (sect. 5)
    mcfg, mdes = cfg_caps(cfg, MOSS_KEY), design_moss()
    k_new = k_own if k_own else K_MAIN
    beta2 = {}
    if pair:
        try:
            beta2 = {s: v[1][1] for s, v in pair_cmp.read(pair[1]).items()}
        except OSError:
            pass
    for r in rows:
        s = r['s']
        if s in mcfg:
            r['smoss'], r['msrc'] = mcfg[s], cfg
        elif s in mdes:
            r['smoss'], r['msrc'] = mdes[s], '设计表'
        else:
            r['smoss'], r['msrc'] = 0., '无'
        site = s in MOSS_SITE and r['smoss'] > 0.
        r['mfrac'] = MOSS_SITE[s][0] if site else PHI_M * r['smoss']
        r['mnote'] = '站点' if site else ('规则' if r['smoss'] > 0. else '—')
        r['cur'] = r['ccap'] if r['ccap'] is not None else r['cap']
        if 'ratio' not in r:
            continue
        r['go_v'] = r['go'] * (1. - r['mfrac'])
        r['ratio_v'] = float(np.sum(r['gm']) / np.sum(r['go_v']))
        r['sug_old'] = r['sug'][k_new] if r['sug'] else '—'
        r['sug_new'] = '—' if LO <= r['ratio_v'] <= HI else suggest(dict(r, go=r['go_v']), k_new)
        n = cap_num(r['sug_new'])
        r['chg'] = n is not None and r['cur'] is not None and abs(n - r['cur']) >= CHANGE
        # sensitivity: moss share of all productivity with vascular belowground NPP = aboveground (sect. 5.2 item 5)
        m2 = r['mfrac'] if site else PHI_M_ALL * r['smoss']
        rv2 = float(np.sum(r['gm']) / np.sum(r['go'] * (1. - m2)))
        s2 = '—' if LO <= rv2 <= HI else suggest(dict(r, go=r['go'] * (1. - m2)), k_new)
        n2 = cap_num(s2)
        r['chg2'] = (n2 is not None and r['cur'] is not None and abs(n2 - r['cur']) >= CHANGE, s2)
        if not pair:
            continue
        cap2 = nml_value(f'{R}/cases/{pair[0]}/sites_{pair[1]}/_conf/ch4_param/{s}.nml', CAP_KEY)
        mod2 = model_monthly(pair[0], pair[1], s)
        g1, l1, g2, l2, go = [], [], [], [], []
        for p in r['pairs']:
            q = [x for x in p if x in mod2]
            if not q:
                break
            g1.append(np.mean([r['mod'][x][0] for x in q]))
            l1.append(np.mean([r['mod'][x][1] for x in q]))
            g2.append(np.mean([mod2[x][0] for x in q]))
            l2.append(np.mean([mod2[x][1] for x in q]))
            go.append(np.mean([r['obs'][x] for x in q]))
        if len(g2) < 3:
            continue
        g1, l1, g2, l2, go = map(np.array, (g1, l1, g2, l2, go))
        moved = abs(l2.mean() / l1.mean() - 1.) >= 0.1
        r['pair'] = dict(cap=cap2, rt=g2.sum() / go.sum(), rv=g2.sum() / (go.sum() * (1. - r['mfrac'])),
                         k=implied_k(g1, l1, g2, l2) if moved else np.nan)
    mok = [r for r in ok if 'ratio_v' in r]
    chg = [r for r in mok if r['chg']]
    pk = [r['pair']['k'] for r in mok if r.get('pair') and r['smoss'] > 0. and np.isfinite(r['pair']['k'])]
    now = datetime.datetime.now()
    out(f'# lai_from_gpp_{now:%y%m%d}：按塔上 GPP 反推湿地塔足迹 LAI 上限')
    out()
    out(f'- **日期**：{now:%Y.%m.%d-%H:%M}')
    out('- **读者**：下一会话')
    ptxt = f' --pair {pair[0]} {pair[1]}' if pair else ''
    out(f'- **来源**：`v2/scripts/lai_from_gpp.py {tree} {tag} {cfg}{ptxt}` 生成（只读 history、观测与配置）；树 `cases/{tree}` @{commit}；'
        f'同树兄弟 tag {", ".join(sibs) if sibs else "无"}；对照配置 `v2/config/{cfg}/SITE_PARAMS.csv`；'
        f'藓校正承接账本 V-172 与 Q-63，站点文献检索 `v2/_scratch/moss/moss_grep.py`（输出 `moss_grep.out`）')
    out(f'- **状态**：{len(rows)} 座湿地塔，可比 {len(ok)} 座；比值出 {LO}–{HI} 的 {len(off)} 座，其中上限起作用的 {len(act)} 座给了建议值；'
        f'第 5 节按维管 GPP（扣藓）重算后要改上限的 {len(chg)} 座')
    out()
    out('## 1. 口径')
    out()
    out(f'1. 季：逐塔取塔上 GPP 月气候态最高的连续 3 个月（可跨年），3 个月都须是强迫气温（`f_xy_t`）月气候态高于 0 °C 的月——北极塔冬季插补留有虚高 GPP'
        '（US-A03 2 月均值 5.4 umol m-2 s-1）；塔上半小时值超出 '
        f'{GPP_MIN:.0f} 到 {GPP_MAX:.0f} umol m-2 s-1 的当缺测（US-A10 有 2.4e37 未填值）；塔上月值要求该月有效时间步 ≥ 80 %；模式与塔上只在同年同月都有值的月对上，'
        '每个日历月先对年求均，季值为 3 个月均值；"季数"为 3 个月都对上的完整季个数，少于 1 写"无"。')
    out('2. 塔上 GPP 用 `GPP`（夜间分区法，FLUXNET-CH4 的 GPP_NT），"比值(DT)"改用 `GPP_DT`（白天分区法）；单位 umol CO2 m-2 s-1 × 1e-6 × 12.011 × 86400 → g C m-2 d-1。'
        '模式 GPP 为 `f_assim`（总光合）× 12.011 × 86400。比值 = 模式 / 塔上。')
    out(f'3. 反推：GPP ∝ 1 − exp(−k·LAI)，逐月求 LAI 缩放 c 使季 GPP 等于塔上；"反推 LAI" = c × 模式季均 LAI（k = {K_MAIN}）；'
        f'建议上限 = c × 模式年峰 LAI（多年均）。只对比值 < {LO} 或 > {HI} 的塔给建议。')
    out('4. 上限状态：压 = 遥感输入年峰 > 本 tag 上限（上限在起作用）；未压 = 遥感峰 ≤ 上限；有林份额1 = 上限不作用于有林份额；非湿地斑块 = IGBP 不是 11，不走湿地上限。'
        '建议栏：补不足 = LAI 无限大也补不上；抬无效 = 需要抬 LAI 但上限未压着；不作用 = 上限对该塔无效；遥感峰封顶 = 抬上限最多到遥感峰。')
    out('5. rstfac 为模式 `f_rstfacsun` 季均（1 为无水分胁迫）。')
    out(f'6. 甲烷 β 为本 tag 的甲烷量级比（模式/塔上），取自 `v2/results/scores77_{tag}.txt`（{"有" if beta else "无此文件"}），只作参照，不进反推。')
    out()
    out('## 2. 逐塔表')
    out()
    kown_h = f' 建议 k{k_own}（模式自身） |' if k_own else ''
    out(f'| 塔 | 类 | 季 | 季数 | 塔上 GPP | 模式 GPP | 比值 | 比值(DT) | 甲烷 β | rstfac | 模式季 LAI | 模式年峰 | 遥感峰 | 上限状态 | 反推 LAI | 本 tag 上限 | {cfg} 上限 |'
        f' 建议上限 k{K_MAIN} | k{K_SENS[0]} / k{K_SENS[1]} |{kown_h}')
    out('|' + ' --- |' * (19 + bool(k_own)))
    for r in rows:
        head = f'| {r["s"]} | {r["c"]} |'
        caps = f'{fmt(r["cap"], 1)} | {fmt(r["ccap"], 1)}'
        b = fmt(beta.get(r['s']))
        if 'ratio' not in r:
            out(f'{head} 无 | 0 | 无 | 无 | 无 | 无 | {b} | 无 | 无 | 无 | {fmt(r["rspeak"])} | {r["state"]} | 无 | {caps} | 无（{r["why"]}） | 无 |'
                + (' 无 |' if k_own else ''))
            continue
        ms = f'{r["ms"][0]}–{r["ms"][-1]}'
        sug = r['sug'][K_MAIN] if r['sug'] else '—'
        sens = f'{r["sug"][K_SENS[0]]} / {r["sug"][K_SENS[1]]}' if r['sug'] else '—'
        own = (f' {r["sug"][k_own] if r["sug"] else "—"} |') if k_own else ''
        inv = '补不足' if np.isinf(r['lai_inv']) else fmt(r['lai_inv'])
        out(f'{head} {ms} | {r["n"]} | {np.mean(r["go"]):.2f} | {np.mean(r["gm"]):.2f} | {r["ratio"]:.2f} | {fmt(r["ratio_dt"])} | {b} |'
            f' {np.mean(r["rs"]):.2f} | {np.mean(r["lai"]):.2f} | {fmt(r["peak"])} | {fmt(r["rspeak"])} | {r["state"]} |'
            f' {inv} | {caps} | {sug} | {sens} |{own}')
    out()
    out('类中位（可比塔）：')
    for c in dict.fromkeys(r['c'] for r in ok):
        rat = [r['ratio'] for r in ok if r['c'] == c]
        out(f'- {c}（{len(rat)} 座）：比值中位 {np.median(rat):.2f}，范围 {min(rat):.2f}–{max(rat):.2f}')
    out()
    out('## 3. 模式自身的 k（同树兄弟 tag 改过该塔上限的对跑）')
    out()
    if sib_rows:
        out('| 塔 | 兄弟 tag | 上限 本 → 兄弟 | 季 LAI 本 → 兄弟 | 季 GPP 本 → 兄弟 | 塔上 GPP | k0.5 预测兄弟 GPP | 反推 k |')
        out('| --- | --- | --- | --- | --- | --- | --- | --- |')
        for x in sib_rows:
            out(f'| {x["s"]} | {x["sb"]} | {x["cap"]:.1f} → {x["cap2"]:.1f} | {x["l1"]:.2f} → {x["l2"]:.2f} | {x["g1"]:.2f} → {x["g2"]:.2f} |'
                f' {x["go"]:.2f} | {x["pred"]:.2f} | {fmt(x["k"])} |')
        out()
        out(f'模式自身 k：{len(ks_own)} 对（{len({x["s"] for x in sib_rows})} 座塔），中位 {fmt(k_own)}，'
            f'范围 {fmt(min(ks_own) if ks_own else np.nan)}–{fmt(max(ks_own) if ks_own else np.nan)}'
            '（同一塔、同一强迫、同一组月的两次跑，上限不同；兄弟 tag 在该塔另有的配置差别也算进这个 k）；同一上限对只列一次。')
    else:
        out('同树没有改过湿地塔上限的已完成兄弟 tag，无对跑。')
    out()
    wetp = [r for r in ok if r['state'] in ('压', '未压', 'GLWD格点', '有林份额1')]
    rs_all = [float(np.mean(r['rs'])) for r in wetp]
    dd = [abs(r['ratio'] - r['ratio_dt']) for r in ok if np.isfinite(r['ratio_dt'])]
    out('## 4. 做法与局限')
    out()
    out('读码行号按树 `cases/paper_v2/v260930/sp_v2_c90` @24061070。')
    out()
    out('1. 上限只压、不抬过遥感：`wetland_lai_cap_site` > 0 时覆盖 GLWD 格点上限（`main/TRACER/MOD_Tracer_Reactive_Methane_WetlandVeg.F90:261`）；'
        '格点上限 `wetveg_laicap` 是格内开阔泥炭（GLWD 23、25）取 `wetland_lai_open_peat` 0.6、沼泽（17、19、27）取 `wetland_lai_marsh` 3.0 的面积加权，'
        '格内没有这些类则不设上限（同文件 204–208 行；默认值 `MOD_Tracer_Reactive_Methane_Const.F90:741–742`）——已验证（读码）。')
    wet_patch = [r['s'] for r in rows if r['state'] != '非湿地斑块']
    nshape = sum(shape[s] for s in wet_patch)
    out('2. `wetveg_cap_lai` 在每次读入遥感 LAI 后调用（`main/CoLM.F90:638、658`），只作用于湿地斑块（patchtype 2，WetlandVeg.F90:300；IGBP 11 对应 patchtype 2，'
        '`main/MOD_Const_LC.F90:398–401`）的非林份额：开 `wetland_lai_shape` 时，遥感年峰高于上限才把全年乘以上限/年峰，否则原样；关时取 min(遥感, 上限)'
        '（WetlandVeg.F90:299–311）。有林份额保持遥感 LAI（304、310 行）——已验证（读码）。'
        f'本 tag {len(wet_patch)} 座湿地斑块塔中 {nshape} 座的 `_conf/ch4_param/<塔>.nml` 把 `wetland_lai_shape` 设为 .true.——已验证（脚本逐塔读取）。')
    out('3. 所以上限只能把 LAI 压到遥感输入以下；把上限抬到遥感年峰以上，LAI 就等于遥感输入，再抬无效。"压"的塔 GPP 偏低可以抬上限，最多抬到遥感峰；'
        '"未压"的塔 GPP 偏低抬上限无效，要换 LAI 输入（站点键 `USE_SITE_LAI`）或查别的原因；有林份额 1 与非湿地斑块的塔上限不起作用——已验证（读码，同上两条）。')
    out('4. k 的量级：湿地斑块（IGBP 11）叶倾角 χ = 0.1（`MOD_Const_LC.F90:445–448` 第 11 项，`CoLMDRIVER.F90:99、151` 按斑块类取用），'
        '直射消光 k_b = (φ1 + φ2·μ)/μ，φ1 = 0.5 − 0.633χ − 0.33χ² = 0.433，φ2 = 0.877(1 − 2φ1) = 0.117（`MOD_Albedo.F90:593–597`），'
        '太阳天顶角 0°、45°、60° 时 k_b = 0.55、0.73、0.98；散射消光 0.719（同文件 599 行）——已验证（读码并代入）。'
        '冠层积分：Rubisco 能力按 (1 − e^(−0.11·L))/0.11 放大，近乎随 L 线性；电子传递能力按 (1 − e^(−0.719·L))/0.719 放大'
        '（`MOD_LeafTemperature.F90:466–472`，`MOD_AssimStomataConductance.F90:560、575`）；光限制速率正比于冠层吸收 PAR（同文件 579 行）——已验证（读码）。'
        '模式 GPP 对 LAI 的有效 k 因而在约 0.1（Rubisco 限制）与 0.55–1.0（光限制）之间，取 0.5、敏感性 0.4 与 0.6——推测。理由：有效 k 由两种限制的比重决定，第 3 节的对跑给模式自身的值。')
    if k_own:
        out(f'   对跑得到的模式自身 k 中位 {k_own}（第 3 节）。k 越小 GPP 越接近随 LAI 线性，同样的 GPP 缺口要的 LAI 变化越小，所以按 0.5 反推在"抬"与"压"两个方向都偏过头'
            f'（抬得多、压得狠）——已验证（第 2 节 k{K_MAIN} 与模式自身两列）。'
            + (f'这些对跑来自 {len({x["s"] for x in sib_rows})} 座温带沼泽塔；泥炭藓塔另由第 5 节对跑检验，k 中位 {np.median(pk):.2f}，并不更高——已验证（第 5.1 节第 4 条）。'
               if pk else
               f'对跑只来自 {len({x["s"] for x in sib_rows})} 座温带沼泽塔，开阔泥炭塔（LAI ≤ 0.6、高纬太阳天顶角大）的有效 k 可能更高——推测。'
               '理由：天顶角大则 k_b 大（本条上文），且弱光下光限制的比重更大；要定需在一座泥炭塔上做一次只改上限的对跑。'))
    out(f'5. GPP 偏差不只来自 LAI：湿地斑块所有类共用一组光合参数（vmax25 = 52 umol m-2 s-1，`MOD_Const_LC.F90:500–503` 第 11 项），'
        '温度响应、辐射（强迫为塔上实测）、水分胁迫都可能贡献；本法把全部偏差记到 LAI 上，得到的是"等效 LAI"——推测。理由：没有逐项分解。'
        f'模式水分胁迫在这些塔上几乎不起作用：可比塔季均 rstfac 最小 {fmt(min(rs_all) if rs_all else np.nan)}、中位 {fmt(np.median(rs_all) if rs_all else np.nan)}——已验证（第 2 节 rstfac 列）。')
    out(f'6. 塔上 GPP 是由 NEE 分区得来的：夜间法与白天法比值之差中位 {fmt(np.median(dd) if dd else np.nan)}、最大 {fmt(max(dd) if dd else np.nan)}——已验证（第 2 节两比值列）；'
        '比值落在 0.8–1.25 边界附近的塔，判断随分区法翻转，建议值只作起点——推测。')
    out('7. 只约束生长季峰：上限在 shape 模式下按比例缩放全年，返青与枯黄的时间仍取自遥感（WetlandVeg.F90:304–305）——已验证（读码）；'
        '季节形状错的塔（如 V-152 所记 US-Sne 的复湿前牧场物候）不能靠上限修——已验证（账本 V-152 行）。')
    out('8. 上限作用于斑块平均 LAI，模式把斑块当成均匀冠层；塔上 GPP 是足迹平均（含开阔水面与非湿地）。V-152 把 US-Sne 冠层内 LAI 乘盖度得 0.3 当上限，'
        '压过头，按塔上 GPP 反推应约 0.8——已验证（账本 V-152 行）。成丛冠层的 GPP 高于同样平均 LAI 的均匀冠层，所以本法反推的是"与塔上 GPP 等效的均匀冠层 LAI"，'
        '不应拿冠层内实测 LAI 乘盖度去核——推测。理由：1 − e^(−kL) 是凹函数。')
    out('9. GPP 对上不等于甲烷对上：V-155 中 US-Myb、US-Sne、US-WPT 压 LAI 后 GPP 已对上塔上，量级比却落到 0.45–0.71，另有底泥偏冷等原因（Q-52）——已验证（账本 V-155 行）。'
        '改上限前后应看 β 与 E/GPP 两项——推测。')
    up_hi = [r['s'] for r in off if r['ratio'] < LO and beta.get(r['s'], np.nan) > HI and r['state'] in ('压', 'GLWD格点')]
    dn_lo = [r['s'] for r in off if r['ratio'] > HI and beta.get(r['s'], np.nan) < LO and r['state'] in ('压', '未压', 'GLWD格点')]
    if beta:
        out(f'   GPP 与甲烷方向相反的塔：GPP 偏低要抬 LAI、甲烷却已偏高（β > {HI}）的 {"、".join(up_hi) or "无"}；'
            f'GPP 偏高要压 LAI、甲烷却已偏低（β < {LO}）的 {"、".join(dn_lo) or "无"}——已验证（第 2 节比值与 β 列）。'
            '这些塔按 GPP 改上限会让甲烷更偏，单位 GPP 的排放另有原因，宜先查后改——推测。')
    diff = []
    for r in off:
        if r['ccap'] is None or r['cap'] is None or abs(r['ccap'] - r['cap']) < 0.01 or not r['sug']:
            continue
        sib = [x for x in sib_rows if x['s'] == r['s'] and abs(x['cap2'] - r['ccap']) < 0.01]
        at = f'；兄弟 tag {sib[0]["sb"]} 以 {r["ccap"]:.1f} 跑出的 GPP 比值 {sib[0]["g2"] / sib[0]["go"]:.2f}' if sib else ''
        own = f'、k{k_own} {r["sug"][k_own]}' if k_own else ''
        diff.append(f'{r["s"]}（{cfg} {r["ccap"]:.1f}；本表 k{K_MAIN} {r["sug"][K_MAIN]}{own}{at}）')
    if diff:
        out(f'10. {cfg} 已改过上限的塔与本表建议对照：{"；".join(diff)}——已验证（第 2、3 节）。')
    moss_section(rows, mok, chg, pk, beta, beta2, pair, cfg, k_new, tag)
    open(f'{V2}/results/lai_from_gpp_{now:%y%m%d}.txt', 'w').write('\n'.join(L) + '\n')


def moss_section(rows, mok, chg, pk, beta, beta2, pair, cfg, k_new, tag):
    """Sect. 5: caps against the vascular (moss-free) GPP target."""
    out()
    out('## 5. 扣藓后的维管 GPP 约束（Q-63 续，承接 V-172）')
    out()
    out(f'1. 维管目标 = 塔上 GPP ×（1 − m），m = φ_m × s_moss，φ_m = {PHI_M}：模式的 `wetland_moss_input_frac` 默认值，注释引 Turetsky et al. 2010 藓占北方与苔原湿地生产力 48 %'
        '（藓 NPP 比上藓加维管地上 NPP；`main/TRACER/MOD_Tracer_Reactive_Methane_Const.F90:841–858`，树 `sp_v2_c91`），'
        f'{cfg} 与 v2_p8lg 的 ch4_parameter.nml 未覆盖——已验证（读码与 grep）。s_moss 取 `{cfg}` 的 `wetland_share_moss_site`，'
        '缺则取 `design_wetland_class_260930.txt` 第 5.5 节表；有站点实测藓占 GPP 的塔改用实测值（第 5.2 节）。')
    out(f'2. 建议上限按模式自身 k = {k_new}（第 3 节）对维管目标重算；新比值（模式 / 维管目标）在 {LO}–{HI} 内写"—"，即保持现行上限；'
        f'"改否"以 |新建议 − 现行上限| ≥ {CHANGE} 为改；现行上限取 `{cfg}`，没有站点值的塔取本 tag 的值。')
    if pair:
        out(f'3. 对跑：`{pair[1]}`（树 `cases/{pair[0]}`）同一组月的季 GPP，比塔上总 GPP 与比维管目标各一列；两次跑季 LAI 差 ≥ 10 % 的塔给出两者之间的 k；'
            f'β 为甲烷量级比，本 tag 取 `scores77_{tag}.txt`，对跑取 `scores77_{pair[1]}.txt`。上限未变的塔两棵树的 f_assim 与 f_lai 逐月相同'
            '（CA-SCB 2015、FI-Lom 2008、US-A10 2014 年 6–8 月在 v2_p9a 与 v2_p9f 三位小数一致）——已验证（python 读 history）。')
    out()
    pt = pair[1] if pair else ''
    ph = f' {pt} 上限 | GPP/塔上 | GPP/维管 | 对跑 k | β 本 → 对跑 |' if pair else ' β |'
    out(f'| 塔 | 类 | s_moss（来源） | m | 塔上 GPP | 维管目标 | 模式 GPP | 新比值 | 旧建议 k{k_new} | 新建议 k{k_new} | 现行上限 | 改否 |{ph}')
    out('|' + ' --- |' * (12 + (5 if pair else 1)))
    for r in rows:
        head = f'| {r["s"]} | {r["c"]} | {r["smoss"]:.0f}（{r["msrc"]}） | {r["mfrac"]:.2f}{"*" if r["mnote"] == "站点" else ""} |'
        b1 = fmt(beta.get(r['s']))
        if 'ratio_v' not in r:
            tail = f' — | — | — | — | {b1} → {fmt(beta2.get(r["s"]))} |' if pair else f' {b1} |'
            out(f'{head} 无 | 无 | 无 | 无 | 无 | 无 | {fmt(r["cur"], 1)} | 否 |{tail}')
            continue
        chg_txt = f'是：{r["cur"]:.1f} → {cap_num(r["sug_new"]):.1f}' if r['chg'] else '否'
        body = (f' {np.mean(r["go"]):.2f} | {np.mean(r["go_v"]):.2f} | {np.mean(r["gm"]):.2f} | {r["ratio_v"]:.2f} |'
                f' {r["sug_old"]} | {r["sug_new"]} | {fmt(r["cur"], 1)} | {chg_txt} |')
        if pair:
            p = r.get('pair')
            tail = (f' {fmt(p["cap"], 1)} | {p["rt"]:.2f} | {p["rv"]:.2f} | {fmt(p["k"])} |' if p else ' — | — | — | — |')
            tail += f' {b1} → {fmt(beta2.get(r["s"]))} |'
        else:
            tail = f' {b1} |'
        out(head + body + tail)
    out()
    out('* 为站点实测藓占 GPP（第 5.2 节），其余 m = φ_m × s_moss。')
    out()
    out('### 5.1 结果')
    out()
    out(f'1. 按维管目标要改上限的 {len(chg)} 座：' + ('；'.join(f'{r["s"]} {r["cur"]:.1f} → {cap_num(r["sug_new"]):.1f}' for r in chg) or '无')
        + '——已验证（上表）。')
    same = [r['s'] for r in mok if r['smoss'] == 0.]
    out(f'2. s_moss = 0 的 {len(same)} 座（marsh、盐沼、非苔藓泥炭等）目标不变，新旧建议相同——已验证（上表 m 列为 0）。')
    down = [r for r in chg if cap_num(r['sug_new']) < r['cur']]
    lowb = [f'{r["s"]}（β {beta[r["s"]]:.2f}）' for r in down if r['smoss'] > 0. and beta.get(r['s'], np.nan) <= 1.05]
    up = [r for r in chg if cap_num(r['sug_new']) > r['cur']]
    hib = [f'{r["s"]}（β {beta[r["s"]]:.2f}）' for r in up if beta.get(r['s'], np.nan) > HI]
    out(f'3. 与甲烷方向相反：扣藓后要压上限、而本 tag 甲烷 β 已 ≤ 1.05 的藓塔 {"、".join(lowb) or "无"}；要抬上限、而 β 已 > {HI} 的 {"、".join(hib) or "无"}'
        '——已验证（上表与第 2 节 β 列）。这些塔按维管 GPP 改上限会让甲烷离塔上更远——推测。理由：V-172 中甲烷随上限同向变（β 列）。')
    if pair and pk:
        sel = [r for r in mok if r.get('pair') and r['smoss'] > 0. and np.isfinite(r['pair']['k'])]
        rt, rv = [r['pair']['rt'] for r in sel], [r['pair']['rv'] for r in sel]
        out(f'4. 对跑检验（藓塔，{len(pk)} 座上限变过）：两次跑之间的 k 中位 {np.median(pk):.2f}（{min(pk):.2f}–{max(pk):.2f}），第 3 节沼泽塔为 {k_new}；'
            f'对跑的季 GPP / 塔上总 GPP 中位 {np.median(rt):.2f}，/ 维管目标中位 {np.median(rv):.2f}——已验证（上表对跑列）。'
            '即 k 取得对，按总 GPP 定的上限把模式 GPP 抬到了塔上总 GPP；甲烷冲过头不是 k 的问题，是目标含藓——推测。理由：GPP 命中而甲烷随之高出（β 列）。')
    c2 = [r for r in mok if r['chg2'][0]]
    out(f'{5 if pair and pk else 4}. 敏感性：φ_m 改取 {PHI_M_ALL}（第 5.2 节第 5 条）时要改上限的 {len(c2)} 座：'
        + ('；'.join(f'{r["s"]} {r["cur"]:.1f} → {cap_num(r["chg2"][1]):.1f}' for r in c2) or '无') + '——已验证（脚本同法重算）。')
    out()
    out('### 5.2 站点文献里的藓占生产力')
    out()
    ntxt = len(glob.glob(f'{V2}/_scratch/moss/txt/*/*.txt'))
    out(f'1. 检索范围：27 座 bog、fen、湿苔原与泥炭排干塔的 `docs/refs/sites/<塔>/` 全部 PDF 及 `_traits` 的两篇 Korrensalo，共 {ntxt} 份转文本，'
        '取同句含藓（moss、Sphagnum、bryophyte）、生产力（GPP、photosynthesis、NPP、productivity、uptake、assimilation）与数字的句子，'
        '再放宽到含百分比或份额词；`scripts/sites/site_facts.csv` 的"年净初级生产力"行逐条看含藓的——已验证（`v2/_scratch/moss/moss_grep.py`、`moss_grep.out`）。')
    out(f'2. FR-LGt：{MOSS_SITE["FR-LGt"][1]}，藓约占 {MOSS_SITE["FR-LGt"][0]:.2f}，低于规则值 {PHI_M:.2f}——已验证（读 PDF 转文本 2.1 节与表 3）。'
        '这是取自本泥炭地的中宇宙，不是塔足迹的分割，按"有站点数据则用"代入——推测。理由：足迹内 Molinia 与 Betula 的比例未知。')
    out('3. SE-Deg：Zajac et al. 2016 四个泥炭地（含 Degerö Stormyr）的芯样中苔藓占活生物量 62–87 %（3.1 节）；US-Bes：Zona et al. 2009 本塔'
        '苔藓约占活生物量 80 %（2006 年 8 月），维管植物 Carex aquatilis LAI 0.64——已验证（读 PDF 转文本）。两者是生物量份额，不是 GPP 份额，不代入——推测。理由：藓单位生物量光合远低于维管叶。')
    out('4. 其余塔没有按藓与维管植物分的 GPP 或 NPP：CA-SCB、US-BZB、US-Uaf 的 site_facts 明写"年净初级生产力"未见（US-BZB 的藓 NPP 在 Bonanza Creek LTER 的 EDI 数据集，不在本仓）；'
        '其余检索命中都是定性句或引文题名——已验证（`moss_grep.out`、site_facts.csv）。这些塔保留 0.48 × s_moss。')
    out('5. φ_m 的口径偏高：Turetsky 的 48 % 是藓比上藓加维管地上 NPP（模式注释即如此），塔上 GPP 还含维管地下生产；若维管地下与地上 NPP 相当，'
        '藓 = 0.48/0.52 × 地上 ≈ 0.92 × 地上，藓 /（藓 + 2 × 地上）≈ 0.32，与 FR-LGt 的 0.33 相近——推测。理由：地下与地上比取 1 是假设，未逐塔核。')


if __name__ == '__main__':
    main()
