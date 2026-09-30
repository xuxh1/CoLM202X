#!/usr/bin/env python
"""Wetland CH4 of a V2 global case by GLWD v2 class (Lehner et al. 2025), mean
over the years given, against per-area literature values by wetland type.

Emission is taken as in global_eval.py (ch4_budget_core land-area weights):
  wetland tile  f_methane_surf_flux_wetland x landarea
  floodplain    f_methane_surf_flux_soil x landarea in cells whose soil tile
                floods in the year (annual max of f_methane_soil_finundated
                >= 1e-4); the rest of the soil tile is upland and left out.
Class areas: v2/data/glwd33_2deg.nc (mk_glwd33.py, km2 per 2-degree cell).
In each cell the wetland-tile emission goes to classes 16-19 and 22-27 (the
classes the tile was built from) and the floodplain emission to classes 8-15
(riverine and lacustrine; they reach the model only through routed
inundation of the soil tile), each in proportion to the class areas in that
cell. A cell with emission but none of those classes counts as unassigned.
Per class and per class group: emission E (Tg CH4/yr), GLWD area A (Mkm2),
E/A (g CH4 m-2 yr-1), then the same by latitude band, a check against
year_budget and the global_eval log, and the literature references.

Writes v2/results/class_eval_<case name>_<y0>-<y1>.log and refuses to
overwrite it. Runs on a compute node, never on a login node.
Usage: class_eval.py <case under cases/paper_v2> <y0> <y1>
  e.g. class_eval.py v260929/g2_o6 2009 2012"""
import os
import re
import socket
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2budget as VB                                                # noqa: E402

B = VB.B
V2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GLWD33 = f'{V2}/data/glwd33_2deg.nc'
LEGEND = '/share/home/dq076/mode/Methane/data/GLWD_v2/area_by_class_pct/GLWD_Legend_v2_0.csv'

TILE = (16, 17, 18, 19, 22, 23, 24, 25, 26, 27)
FLOOD = tuple(range(8, 16))
GROUPS = (('boreal peat', (22, 23)), ('temperate peat', (24, 25)), ('tropical peat', (26, 27)),
          ('palustrine regularly flooded', (16, 17)), ('palustrine seasonally saturated', (18, 19)),
          ('riverine + lacustrine', FLOOD))
BANDS = (('90S-30S', -90, -30), ('30S-30N', -30, 30), ('30N-60N', 30, 60), ('60N-90N', 60, 90))

# Literature per class group, g CH4 m-2 yr-1 unless marked. Values and DOIs are
# copied from the surveys named in each line (their DOIs were checked against
# Crossref / DataCite there); g CH4-C is converted with x 1.336.
REFS = {
    'boreal peat': [
        'BAWLD-CH4, Kuhn et al. 2021 Table 3 p.5170, WARM-SEASON site fluxes in mg CH4 m-2 d-1 (not annual), '
        'median [IQR]: permafrost bog 2.32 [0-6.9] n81; bog 24.55 [6.92-57.35] n87; fen 54 [20-107.2] n109; '
        'tundra wetland 65 [34-99.3] n109',
        'Kuhn et al. 2025, annual per class (Suppl. Table S1 / class area of Olefeldt et al. 2021 Table 3), no IQR: '
        'permafrost bog ~1.2; bog ~4.2; fen ~12.5; tundra wetland ~11.5; BAWLD area-weighted mean 7.55',
        'Treat et al. 2018 data set, boreal annual, median [IQR]: bog 5.6 [3.2-13.1] n110; fen 8.4 [4.4-18.7] n227; '
        'measured-only annual medians bog 6.8 n25, fen 17.6 n39; tundra bog 0.9',
        'Treat et al. 2018 (GCB) p.13, biome annual median: tundra 6.2 +- 1.7; boreal 7.2 +- 1.4',
    ],
    'temperate peat': [
        'BAWLD-CH4: none (boreal-Arctic domain only)',
        'Treat et al. 2018 (GCB) p.13, biome annual median: temperate 13.3 +- 5.4; IQR and bog/fen split: none',
    ],
    'tropical peat': [
        'Treat et al. 2018: none (temperate, boreal and Arctic sites only); median and IQR: none published',
        'site values (lit_tropical_mech_260929.txt sec 2): SE Asia intact peat swamp forest, chamber meta-analysis 3.9 '
        "(2.9 g C; Hergoualc'h & Verchot 2014 via Hergoualc'h et al. 2020 p.7212; Griffis et al. 2020 prints 28.6 g C "
        'for the same, unresolved); SE Asia EC MY-MLM 12.7, Kampar 9.1, ID-Pag 0.01-0.23',
        "site values (same survey): Peru Quistococha EC 29.4 [26.7-32.1] (Griffis et al. 2020), intact-forest chamber "
        "30.2 (Hergoualc'h et al. 2020 Table 6), PMFB diffusive 17.6 (Teh et al. 2017); Congo peat: no in situ value",
    ],
    'palustrine regularly flooded': [
        'BAWLD-CH4 marsh (boreal-Arctic), Kuhn et al. 2021 Table 3: WARM-SEASON 106 [70.5-200] mg CH4 m-2 d-1 n33; '
        'Kuhn et al. 2025 marsh annual ~33',
        'Treat et al. 2018 data set, boreal marsh annual median 23.5 n21 (IQR: none in the survey)',
        'tropical permanent papyrus swamp, one tower: BW-Gum 115.9 (2018), 122.1 (2019) (Helfter et al. 2022 '
        'Phil Trans A Table 1; lit_floodplain_260929.txt sec 1.1)',
    ],
    'palustrine seasonally saturated': [
        'none',
    ],
    'riverine + lacustrine': [
        'BAWLD-CH4, Treat et al. 2018: no floodplain class (none); median and IQR: none',
        'tropical seasonally flooded, annual per land area: BR-Npw Pantanal flooded forest 25.7 (19.21 g C; Delwiche '
        'et al. 2021 Table B3(c)); BW-Nxr Okavango floodplain 36.3 (2018), 8.7 (2019 drought) (Helfter et al. 2022 '
        'Phil Trans A Table 1); Congo flooded forest 16.8-33.5 (12.6-25.1 g C; Tathy et al. 1992, abstract)',
    ],
}
REF_ALL = ('tropical wetlands, all types (not annual): compilation of 328 points, median 35 [5-160] mg CH4 m-2 d-1 '
           '(Murguia-Flores et al. 2023, abstract; lit_floodplain_260929.txt sec 1.5)')
DOIS = ('Kuhn 2021 10.5194/essd-13-5151-2021; Kuhn 2025 10.1038/s41558-025-02413-y; '
        'Olefeldt 2021 10.5194/essd-13-5127-2021; Treat 2018 data set 10.1594/PANGAEA.886976; '
        'Treat 2018 GCB 10.1111/gcb.14137; Hergoualc\'h 2020 10.1111/gcb.15354; '
        'Hergoualc\'h & Verchot 2014 10.1007/s11027-013-9511-x (not read); Griffis 2020 10.1016/j.agrformet.2020.108167; '
        'Teh 2017 10.5194/bg-14-3669-2017; Tathy 1992 10.1029/90JD02555; Delwiche 2021 10.5194/essd-13-3607-2021; '
        'Helfter 2022 10.1098/rsta.2021.0148; Murguia-Flores 2023 10.1029/2022GB007601; '
        'MY-MLM, Kampar, ID-Pag: earlier-round notes, no DOI in the survey')


def glwd_on(lat, lon):
    """{class: (nlat, nlon) km2} of GLWD classes 8-19, 22-27 on the model grid (0 outside 84N-56S)."""
    with xr.open_dataset(GLWD33) as g:
        if not np.allclose(g.lon.values, lon):
            sys.exit('glwd33_2deg.nc longitudes differ from the model grid')
        a = {}
        for k in FLOOD + TILE:
            v = g[f'area_class_{k:02d}']
            r = v.reindex(lat=lat, method='nearest', tolerance=0.01).fillna(0).values
            if abs(r.sum() - float(v.sum())) > 1e-6 * max(float(v.sum()), 1.0):
                sys.exit(f'class {k}: area lost when put on the model grid')
            a[k] = r
    return a


def share(stack):
    """Per-cell class shares of a (nclass, nlat, nlon) area stack; 0 where the stack is empty."""
    tot = stack.sum(0)
    return np.divide(stack, tot, out=np.zeros_like(stack), where=tot > 0), tot > 0


def global_eval_log(label, name, y0, y1):
    """(path, (E_wetland, E_floodplain) or a reason string) from the global_eval log of the same window."""
    f = f'{V2}/results/global_eval_{name}_{y0}-{y1}.log'
    if not os.path.exists(f):
        return f, 'no log'
    txt = open(f).read()
    if not txt.startswith(f'# {label},'):
        return f, 'log header names another case or window'
    et = re.search(r'^E_wetland\s+mean\s+(-?[\d.]+)', txt, re.M)
    fp = re.search(r'^E_floodplain .*? mean\s+(-?[\d.]+)', txt, re.M)
    if not (et and fp):
        return f, 'E_wetland or E_floodplain line missing'
    return f, (float(et.group(1)), float(fp.group(1)))


def fmt(E, A):
    return f'{E:7.2f} {A:7.3f} ' + (f'{E / A:7.1f}' if A > 0 else f'{"-":>7s}')


def main():
    if socket.gethostname().startswith('mgt'):
        sys.exit('run on a compute node, not on a login node')
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sub, name = os.path.split(sys.argv[1].strip('/'))
    ver, y0, y1 = f'paper_v2/{sub}', int(sys.argv[2]), int(sys.argv[3])
    out = f'{V2}/results/class_eval_{name}_{y0}-{y1}.log'
    if os.path.exists(out):
        sys.exit(f'{out} exists; remove it first')
    case = B.Case(ver, name)
    years = [y for y in B.available_years(case) if y0 <= y <= y1]
    if not years:
        sys.exit(f'{case.dir}: no tracer history in {y0}-{y1}')
    areas = VB.load_areas(case, years[0])
    A_land = areas[0]

    acc = {k: [] for k in ('et', 'ef', 'at', 'af', 'wet', 'fp')}
    for y in years:
        bud = B.year_budget(y, case, areas)
        with xr.open_dataset(f'{case.dir}/history/{name}_hist_tracer_{y}.nc') as tr:
            w = xr.DataArray(B.SEC_PER_MONTH[:tr.sizes['time']], dims='time')

            def tg(v):
                return np.nan_to_num(((tr[v].fillna(0) * A_land * w).sum('time') * B.M_CH4).values)

            fsoil = tr['f_methane_soil_finundated'].fillna(0)
            dry = (fsoil.max('time') < 1e-4).values
            es = tg('f_methane_surf_flux_soil')
            acc['et'].append(tg('f_methane_surf_flux_wetland'))
            acc['ef'].append(np.where(dry, 0.0, es))
            acc['at'].append(np.nan_to_num((tr['f_methane_area_wetland'].fillna(0).max('time') * A_land).values) / 1e12)
            acc['af'].append(np.nan_to_num(((fsoil * tr['f_methane_area_soil'].fillna(0)).max('time')
                                            * A_land).values) / 1e12)
            # global_eval.py: floodplain = E_soil of year_budget minus the soil tile of never-flooding cells
            acc['wet'].append(bud['E_wetland'])
            acc['fp'].append(bud['E_soil'] - float(np.where(dry, es, 0.0).sum()))
            lat, lon = tr['lat'].values, tr['lon'].values
    m = {k: np.mean(np.array(v), 0) for k, v in acc.items()}

    leg = pd.read_csv(LEGEND).set_index('GLWD_ID')['Class_name']
    a = glwd_on(lat, lon)
    shT, hasT = share(np.stack([a[k] for k in TILE]))
    shF, hasF = share(np.stack([a[k] for k in FLOOD]))
    E = {k: m['et'] * shT[i] for i, k in enumerate(TILE)}
    E.update({k: m['ef'] * shF[i] for i, k in enumerate(FLOOD)})
    un_t, un_f = np.where(hasT, 0.0, m['et']), np.where(hasF, 0.0, m['ef'])
    lat2 = np.broadcast_to(lat[:, None], m['et'].shape)
    Mk = {k: v / 1e6 for k, v in a.items()}                           # Mkm2 per cell

    label = f'{ver}/{name}'
    L = [f'# {label}, {years[0]}-{years[-1]} ({len(years)} years); wetland-tile CH4 split over GLWD v2 classes '
         f'16-19, 22-27 and floodplain CH4 over classes 8-15, by class-area share per 2-degree cell',
         '# E Tg CH4/yr; A GLWD class area Mkm2 (v2/data/glwd33_2deg.nc); E/A g CH4 m-2 yr-1',
         '== per class',
         f'{"cls":>3s} {"E":>7s} {"A":>7s} {"E/A":>7s}  name']
    for k in FLOOD + TILE:
        L.append(f'{k:3d} {fmt(E[k].sum(), Mk[k].sum())}  {leg[k]}')
    L.append(f'unassigned wetland tile (no class 16-19, 22-27 in the cell) {un_t.sum():7.2f}; '
             f'unassigned floodplain (no class 8-15 in the cell) {un_f.sum():7.2f}')

    tot = m['et'].sum() + m['ef'].sum()
    L += ['== per class group',
          f'{"group":32s} {"classes":>8s} {"E":>7s} {"A":>7s} {"E/A":>7s} {"E share":>7s}']
    for g, cls in GROUPS:
        e, ar = sum(E[k].sum() for k in cls), sum(Mk[k].sum() for k in cls)
        L.append(f'{g:32s} {f"{cls[0]}-{cls[-1]}":>8s} {fmt(e, ar)} {100 * e / tot:6.1f}%')
    L.append(f'{"unassigned tile + floodplain":32s} {"":8s} {un_t.sum() + un_f.sum():7.2f}')
    L.append(f'{"all (tile + floodplain)":32s} {"":8s} {tot:7.2f} '
             f'{sum(Mk[k].sum() for k in TILE + FLOOD):7.3f}')

    L += ['== per class group by latitude band (E, A, E/A)',
          f'{"group":32s} ' + ' | '.join(f'{b:^23s}' for b, *_ in BANDS)]
    for g, cls in GROUPS + (('unassigned', ()),):
        cells = []
        for b, lo, hi in BANDS:
            bm = (lat2 >= lo) & (lat2 < hi)
            if cls:
                e = sum(E[k][bm].sum() for k in cls)
                ar = sum(Mk[k][bm].sum() for k in cls)
            else:
                e, ar = un_t[bm].sum() + un_f[bm].sum(), 0.0
            cells.append(fmt(e, ar))
        L.append(f'{g:32s} ' + ' | '.join(cells))
    L.append(f'{"all tile + floodplain":32s} ' + ' | '.join(
        f'{(m["et"] + m["ef"])[(lat2 >= lo) & (lat2 < hi)].sum():7.2f}{"":16s}' for b, lo, hi in BANDS))

    L.append('== check')
    et_sum = sum(E[k].sum() for k in TILE) + un_t.sum()
    ef_sum = sum(E[k].sum() for k in FLOOD) + un_f.sum()
    L.append(f'wetland tile: classes + unassigned {et_sum:.4f} vs year_budget E_wetland {m["wet"]:.4f} '
             f'(diff {et_sum - m["wet"]:+.2e})')
    L.append(f'floodplain:   classes + unassigned {ef_sum:.4f} vs year_budget E_soil - upland {m["fp"]:.4f} '
             f'(diff {ef_sum - m["fp"]:+.2e})')
    f, ge = global_eval_log(label, name, y0, y1)
    if isinstance(ge, tuple):
        L.append(f'global_eval log {os.path.basename(f)}: E_wetland {ge[0]:.1f}, E_floodplain {ge[1]:.1f}; '
                 f'here {et_sum:.2f}, {ef_sum:.2f} (diff {et_sum - ge[0]:+.2f}, {ef_sum - ge[1]:+.2f}; log rounds to 0.1)')
    else:
        L.append(f'global_eval log {os.path.basename(f)}: {ge}')

    L.append('== areas (Mkm2)')
    gT, gF = sum(Mk[k] for k in TILE), sum(Mk[k] for k in FLOOD)
    L.append(f'model wetland tile (static) {m["at"].sum():.3f}; GLWD 16-19, 22-27 {gT.sum():.3f}, of it in cells '
             f'without a model wetland tile {gT[m["at"] <= 0].sum():.3f}')
    L.append(f'model floodplain, sum of per-cell annual max of soil finundated x soil tile {m["af"].sum():.3f}; '
             f'GLWD 8-15 {gF.sum():.3f}, of it in cells whose soil tile never floods {gF[m["af"] <= 0].sum():.3f}')

    L.append('== literature per class group (g CH4 m-2 yr-1 unless marked; g CH4-C x 1.336)')
    for g, _ in GROUPS:
        L.append(f'{g}:')
        L += [f'  - {r}' for r in REFS[g]]
    L += [f'all tropical: {REF_ALL}', f'DOI: {DOIS}']

    txt = '\n'.join(L)
    print(txt)
    with open(out, 'x') as fo:
        fo.write(txt + '\n')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
