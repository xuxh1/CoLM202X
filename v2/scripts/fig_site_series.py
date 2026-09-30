"""Monthly CH4 at the 23 main towers: observation against two V2 site tags.

Pairs and units are score_series' (monthly means of daily FCH4_f with at least
10 valid days; model f_methane_surf_flux_tot_active), so the lines are the
series the KGE terms are computed on. Three pages of 2 x 4 panels, towers in
class order; page geometry follows ~/.ai_rules/draw_colm.md (177 mm wide,
panel rows of 1 U = 46 mm, Arial 12/10/9 pt, data lines 1.4 pt). Each panel
prints KGE_ln | beta | alpha | r of both tags from their summary files.
Writes site_series_<tagA>_vs_<tagB>_<n>.png and a layout log into the project
figures tree at the case's path (figures/paper_v2/<generation>/<case>/, OPS-5),
or into [out dir] when given.
Usage: fig_site_series.py <case dir> <tagA> <labelA> <tagB> <labelB> [out dir]"""
import os
import re
import sys
from datetime import date

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                      # noqa: E402
import matplotlib.dates as mdates                                    # noqa: E402

sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/sites')
sys.path.insert(0, '/share/home/dq076/mode/Methane/scripts/figures')
import score_series as SS                                            # noqa: E402
import ch4_style as ST                                               # noqa: E402

RES = '/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/results'
CLASS_CODE = {'permafrost': 1, 'bog': 3, 'fen': 4, 'marsh': 5, 'salt_marsh': 6, 'trop_swamp': 7}
ORDER = ['permafrost', 'bog', 'fen', 'marsh', 'trop_swamp', 'salt_marsh']
MM = 1 / 25.4
PT = 25.4 / 72                                   # mm per point
W_MM, U = 177.0, 46.0
M_LEFT = M_TOP = M_BOT = ROW_GAP = 0.6 * 9 * PT  # 1.9 mm
M_RIGHT, COL_GAP = 3.8, 9 * PT
TITLE, TICKLAB = 10 * PT, 9 * PT
TICK = 3.5 * PT
YLAB_W = 6.8                                     # widest y tick label, mm
LEFT_DECOR = 10 * PT + 2 * PT + YLAB_W + 2 * TICK  # y title + gap + labels + pad + tick
INNER_DECOR = YLAB_W + 2 * TICK
PLOT_W = (W_MM - M_LEFT - M_RIGHT - LEFT_DECOR - INNER_DECOR - COL_GAP) / 2
PLOT_H = U - TITLE - 3 * PT - 2 * TICK - TICKLAB


def summary(tag):
    """site -> (class, KGEln, beta, alpha, r) from the tag's summary file."""
    out = {}
    for line in open(f'{RES}/{tag}_summary.txt'):
        p = line.split()
        if len(p) == 6 and re.match(r'^[A-Z]{2}-\w{3}$', p[0]):
            out[p[0]] = (p[1], *map(float, p[2:]))
    return out


def series(d):
    """{(y, m): v} -> (dates, values) on a continuous monthly axis, NaN in gaps."""
    if not d:
        return np.array([]), np.array([])
    y0, y1 = min(k[0] for k in d), max(k[0] for k in d)
    ks = [(y, m) for y in range(y0, y1 + 1) for m in range(1, 13)]
    return (np.array([date(y, m, 15) for y, m in ks]),
            np.array([d.get(k, np.nan) for k in ks], dtype=float))


def thin(n, k=16):
    return max(1, int(np.ceil(n / k)))


def style():
    plt.rcParams.update({
        'font.family': 'Arial', 'font.size': 9, 'axes.titlesize': 10, 'axes.titleweight': 'bold',
        'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
        'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
        'xtick.major.size': 3.5, 'ytick.major.size': 3.5, 'xtick.major.pad': 3.5,
        'ytick.major.pad': 3.5, 'xtick.direction': 'out', 'ytick.direction': 'out',
        'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.labelpad': 2,
        'axes.titlepad': 3, 'text.color': 'black', 'axes.edgecolor': 'black',
        'pdf.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})


def panel(ax, sid, cls, obs, ma, mb, sa, sb, la, lb, first):
    _, colour, marker = ST.class_style(CLASS_CODE[cls])
    to, vo = series(obs)
    ta, va = series(ma)
    tb, vb = series(mb)
    for t, v, ls, lab in ((tb, vb, (0, (3, 1.6)), lb), (ta, va, '-', la)):
        if t.size:
            st = thin(np.isfinite(v).sum())
            ax.plot(t, v, color=colour, lw=1.4, ls=ls, marker=marker, ms=5.2,
                    markevery=st, mec='white', mew=0.4, label=lab, zorder=3)
    if to.size:
        st = thin(np.isfinite(vo).sum())
        ax.plot(to, vo, color='#0b0b0b', lw=1.4, marker='o', ms=4.6, markevery=st,
                label='Observation', zorder=4)
    vals = np.concatenate([x[np.isfinite(x)] for x in (vo, va, vb) if x.size])
    lo, hi = min(0.0, vals.min()), vals.max()
    ax.set_ylim(lo, hi + (1.0 if first else 0.55) * (hi - lo))   # room for the texts
    if lo < 0:
        ax.axhline(0, color='#8a8880', ls='--', lw=0.7, zorder=1)
    ax.grid(axis='y', ls=':', lw=0.5, alpha=0.5, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    tt = [t for t in (to, ta, tb) if t.size]
    t0, t1 = min(t[0] for t in tt), max(t[-1] for t in tt)
    ax.set_xlim(date(t0.year, 1, 1), date(t1.year, 12, 31))
    ny = t1.year - t0.year + 1
    if ny <= 2:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, _: (mdates.num2date(x).strftime('%b.') if mdates.num2date(x).month in (1, 4, 7, 10) else '')))
    else:
        step = 1 if ny <= 6 else 2
        ax.xaxis.set_major_locator(mdates.YearLocator(step))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    txt = [f'{la}: {sa[1]:.2f} | {sa[2]:.2f} | {sa[3]:.2f} | {sa[4]:.2f}' if sa else f'{la}: n/a',
           f'{lb}: {sb[1]:.2f} | {sb[2]:.2f} | {sb[3]:.2f} | {sb[4]:.2f}' if sb else f'{lb}: n/a']
    if first:
        txt.insert(0, r'KGE$_\mathrm{ln}$ | $\beta$ | $\alpha$ | $r$')
    ax.text(0.02, 0.98, '\n'.join(txt).replace('-', '−'), transform=ax.transAxes,
            va='top', ha='left', fontsize=9, linespacing=1.15, zorder=5)
    if first:
        # line style names the tag; colour is the tower's class in every panel
        from matplotlib.lines import Line2D
        hs = [Line2D([], [], color='#0b0b0b', lw=1.4, marker='o', ms=4.6, label='Observation'),
              Line2D([], [], color='#8a8880', lw=1.4, ls='-', label=la),
              Line2D([], [], color='#8a8880', lw=1.4, ls=(0, (3, 1.6)), label=lb)]
        ax.legend(handles=hs, loc='upper left', bbox_to_anchor=(0.0, 0.64), ncol=3,
                  frameon=False, handlelength=2.2, columnspacing=1.2, borderaxespad=0.2)


def decor(fig, ax):
    """(left, right, top, bottom) reach of the panel's texts beyond its plot area, mm."""
    rend = fig.canvas.get_renderer()
    tb, bb = ax.get_tightbbox(rend), ax.get_window_extent(rend)
    px = fig.dpi / 25.4
    return ((bb.x0 - tb.x0) / px, max(0.0, tb.x1 - bb.x1) / px,
            (tb.y1 - bb.y1) / px, (bb.y0 - tb.y0) / px)


def place(fig, axes, h_mm, passes=4):
    """Put the plot areas so the measured texts meet the draw_colm margins: outer
    margins 1.9 mm (right 3.8), column gap 3.2 mm text to text, row gap 1.9 mm,
    each panel 1 U tall from its title top to its tick-label bottom; plot areas
    of a column share left and right edges, of a row share top and bottom."""
    ncol = 1 + max(c for _, c, _ in axes)
    nrow = 1 + max(r for r, _, _ in axes)
    for _ in range(passes):
        fig.canvas.draw()
        d = {(r, c): decor(fig, ax) for r, c, ax in axes}
        L = [max(d[k][0] for k in d if k[1] == c) for c in range(ncol)]
        R = [max(d[k][1] for k in d if k[1] == c) for c in range(ncol)]
        T = [max(d[k][2] for k in d if k[0] == r) for r in range(nrow)]
        B = [max(d[k][3] for k in d if k[0] == r) for r in range(nrow)]
        w = (W_MM - M_LEFT - M_RIGHT - sum(L) - sum(R) - COL_GAP * (ncol - 1)) / ncol
        for r, c, ax in axes:
            x0 = M_LEFT + sum(L[:c + 1]) + sum(R[:c]) + c * (w + COL_GAP)
            top = h_mm - M_TOP - r * (U + ROW_GAP) - T[r]
            h = U - T[r] - B[r]
            ax.set_position([x0 / W_MM, (top - h) / h_mm, w / W_MM, h / h_mm])


def check(fig, axes, h_mm, ip, path):
    """draw_colm section 11: measured margins, canvas bounds and overlaps."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px = fig.dpi / 25.4
    tb = [ax.get_tightbbox(rend) for _, _, ax in axes]
    out = []
    for (r, c, ax), b in zip(axes, tb):
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > fig.bbox.width + 0.5 or b.y1 > fig.bbox.height + 0.5:
            out.append(f'page {ip}: {ax.get_title(loc="left")} leaves the canvas')
    for i in range(len(tb)):
        for j in range(i + 1, len(tb)):
            if tb[i].overlaps(tb[j]):
                out.append(f'page {ip}: panels {i} and {j} overlap')
    gaps = []
    for (r, c, ax), b in zip(axes, tb):
        for (r2, c2, ax2), b2 in zip(axes, tb):
            if r2 == r and c2 == c + 1:
                gaps.append((b2.x0 - b.x1) / px)
    rows = sorted({r for r, _, _ in axes})
    rg = [(min(b.y0 for (r, _, _), b in zip(axes, tb) if r == rr) -
           max(b.y1 for (r, _, _), b in zip(axes, tb) if r == rr + 1)) / px
          for rr in rows[:-1] if any(r == rr + 1 for r, _, _ in axes)]
    wid = sorted({round(ax.get_position().width * W_MM, 3) for _, _, ax in axes})
    out.append(f'page {ip}: {len(axes)} panels; margins left {min(b.x0 for b in tb) / px:.2f} '
               f'right {(fig.bbox.width - max(b.x1 for b in tb)) / px:.2f} '
               f'top {(fig.bbox.height - max(b.y1 for b in tb)) / px:.2f} '
               f'bottom {min(b.y0 for b in tb) / px:.2f} mm; column gap '
               f'{min(gaps):.2f}-{max(gaps):.2f} mm; row gap {min(rg):.2f}-{max(rg):.2f} mm; '
               f'plot widths {wid} mm -> {path}')
    return out


def main():
    case, ta, la, tb, lb = sys.argv[1:6]
    out = sys.argv[6] if len(sys.argv) > 6 else os.path.join(
        '/share/home/dq076/mode/Methane/figures',
        os.path.relpath(os.path.abspath(case), '/share/home/dq076/mode/Methane/cases'))
    os.makedirs(out, exist_ok=True)
    style()
    sa, sb = summary(ta), summary(tb)
    sites = sorted(sa, key=lambda s: (ORDER.index(sa[s][0]), s))
    pages = [sites[i:i + 8] for i in range(0, len(sites), 8)]
    h_mm = M_TOP + 4 * U + 3 * ROW_GAP + M_BOT
    log = [f'canvas {W_MM:.1f} x {h_mm:.1f} mm; plot area {PLOT_W:.2f} x {PLOT_H:.2f} mm']
    letters = 'abcdefghijklmnopqrstuvwxyz'
    k = 0
    for ip, page in enumerate(pages, 1):
        fig = plt.figure(figsize=(W_MM * MM, h_mm * MM))
        axes = []
        for j, sid in enumerate(page):
            r, c = divmod(j, 2)
            ax = fig.add_axes([0.1, 0.1, 0.3, 0.1])
            cls = sa[sid][0]
            name = ST.class_style(CLASS_CODE[cls])[0]
            panel(ax, sid, cls, SS.obs_monthly(sid), SS.mod_monthly(case, ta, sid),
                  SS.mod_monthly(case, tb, sid), sa.get(sid), sb.get(sid), la, lb, j == 0)
            ax.set_title(f'({letters[j]}) {sid}, {name}', loc='left')
            if c == 0:
                ax.set_ylabel(r'$\mathrm{CH_4}$ flux (mg m$^{-2}$ d$^{-1}$)')
            axes.append((r, c, ax))
        place(fig, axes, h_mm)
        p = f'{out}/site_series_{ta}_vs_{tb}_{ip}.png'
        fig.savefig(p, dpi=600)
        log += check(fig, axes, h_mm, ip, p)
        plt.close(fig)
    open(f'{out}/site_series_{ta}_vs_{tb}.log', 'w').write('\n'.join(log) + '\n')
    print('\n'.join(log))


if __name__ == '__main__':
    main()
