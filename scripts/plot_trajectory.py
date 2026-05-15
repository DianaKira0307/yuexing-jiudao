"""月行九道 — (λ, β) 天球轨迹可视化

默认读取 data/traj_ecliptic_date.txt，以四个时间窗口展示月球在黄道坐标中走过的路径带。
窗口宽度从 1 个月到整个节点退行周期，直接呈现"月行九道"的几何含义。

输出: data/obs_traj_ecliptic_date.png
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FRAME = 'j2000' if '--j2000' in sys.argv else 'date'
DATA_FILE = f'data/traj_ecliptic_{FRAME}.txt'
OUT_FILE  = f'data/obs_traj_ecliptic_{FRAME}.png'

WINDOWS = [
    ('1 month  (~1 orbit)',              30),
    ('1 year   (~13 orbits)',           365),
    ('9.3 years  (half nodal cycle)',  3397),
    ('18.6 years (full nodal cycle)',  6793),
]


def main():
    if not os.path.exists(DATA_FILE):
        sys.exit(f'[Error] {DATA_FILE} not found. Run: make trajectory')

    print(f'Loading {DATA_FILE} ...')
    data = np.loadtxt(DATA_FILE, comments='#')
    days, lam, beta = data[:, 0], data[:, 1], data[:, 2]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    frame_label = 'SOFA date ecliptic' if FRAME == 'date' else 'fixed J2000 ecliptic'
    fig.suptitle("Lunar path on the celestial sphere — Nine Paths of the Moon\n"
                 f"({frame_label}, DE440 starting 1900-01-01)", fontsize=12)

    colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd']

    for ax, (label, n_days), color in zip(axes.flat, WINDOWS, colors):
        mask = days <= n_days
        ax.scatter(lam[mask], beta[mask],
                   s=1.5, c=color, alpha=0.5, rasterized=True)
        ax.axhline(0,   color='gray', lw=0.6, ls='--', alpha=0.7)
        ax.axhline( 5.1, color='gray', lw=0.4, ls=':', alpha=0.4)
        ax.axhline(-5.1, color='gray', lw=0.4, ls=':', alpha=0.4)
        ax.set_xlim(0, 360)
        ax.set_ylim(-7.5, 7.5)
        ax.set_xlabel('Ecliptic longitude λ (deg)')
        ax.set_ylabel('Ecliptic latitude β (deg)')
        ax.set_title(label)
        ax.grid(True, alpha=0.2)
        n_pts = mask.sum()
        ax.text(0.02, 0.97, f'N = {n_pts}',
                transform=ax.transAxes, va='top', fontsize=8, color='#555')

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f'Saved: {OUT_FILE}')


if __name__ == '__main__':
    main()
