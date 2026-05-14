"""Lunar Nine Paths -- ecliptic power spectrum visualization

Reads data/peaks_{r,dlambda,beta}.txt produced by `make analyze`,
plots the normalized power spectra for geocentric distance r,
longitude residual delta-lambda, and ecliptic latitude beta.

Output: data/spectrum_ecliptic.png
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

DATA_DIR = 'data'
OUT_FILE = 'data/spectrum_ecliptic.png'
N_LABEL  = 8      # number of peaks to annotate per panel
P_MAX    = 40000  # ignore periods beyond this (removes DC artifact)

# (variable name, panel title, stem color)
PANELS = [
    ('r',       'Geocentric distance r  (km)',         '#1f77b4'),
    ('dlambda', 'Longitude residual  d-lambda  (rad)', '#2ca02c'),
    ('beta',    'Ecliptic latitude  beta  (rad)',       '#d62728'),
]


def load_peaks(varname):
    path = os.path.join(DATA_DIR, f'peaks_{varname}.txt')
    if not os.path.exists(path):
        sys.exit(f'[Error] {path} not found.  Run: make analyze')
    periods, powers = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            if len(cols) >= 5:
                p  = float(cols[3])
                pw = float(cols[4])
                if p < P_MAX:
                    periods.append(p)
                    powers.append(pw)
    return np.array(periods), np.array(powers)


def fmt_period(days):
    if days < 200:
        return f'{days:.1f}d'
    elif days < 3650:
        return f'{days/365.25:.2f}yr'
    else:
        return f'{days/365.25:.1f}yr'


def _xfmt(x, _):
    if x < 365:
        return f'{x:.0f}d'
    elif x < 3650:
        return f'{x/365.25:.1f}yr'
    else:
        return f'{x/365.25:.0f}yr'


def draw_panel(ax, varname, title, color):
    periods, powers = load_peaks(varname)
    if len(periods) == 0:
        return

    rel   = powers / powers.max()          # normalize to dominant peak
    order = np.argsort(periods)
    p_plot = periods[order]
    r_plot = rel[order]

    # stem lines + scatter dots
    for p, r in zip(p_plot, r_plot):
        ax.plot([p, p], [1e-8, r], color=color, lw=0.9, alpha=0.6, zorder=2)
    ax.scatter(p_plot, r_plot, s=18, color=color, zorder=3)

    # annotate top N peaks by power
    for p, r in zip(periods[:N_LABEL], rel[:N_LABEL]):
        ax.annotate(
            fmt_period(p),
            xy=(p, r), xytext=(3, 6),
            textcoords='offset points',
            fontsize=7.5, color='#cc2222', fontweight='bold',
            va='bottom', ha='left', zorder=5,
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(6, P_MAX)
    ax.set_ylim(5e-8, 5)
    ax.set_xlabel('Period', fontsize=9)
    ax.set_ylabel('Relative power  (log)', fontsize=9)
    ax.set_title(title, fontsize=10, pad=4)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_xfmt))
    ax.grid(True, which='major', ls='-',  alpha=0.12)
    ax.grid(True, which='minor', ls='--', alpha=0.06)


def main():
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    fig.subplots_adjust(hspace=0.45)
    fig.suptitle(
        'Lunar ecliptic power spectra  (DE440, 1900-2079,  N=65536,  Hanning window)',
        fontsize=11, y=0.99,
    )

    for ax, (varname, title, color) in zip(axes, PANELS):
        draw_panel(ax, varname, title, color)

    fig.savefig(OUT_FILE, dpi=150, bbox_inches='tight')
    print(f'Saved: {OUT_FILE}')


if __name__ == '__main__':
    main()
