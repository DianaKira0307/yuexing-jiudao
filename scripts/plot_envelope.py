"""月行九道 — β 振幅包络可视化

默认读取 data/envelope_beta_date.txt，绘制月球黄纬 β 的上下包络。
包络幅度的缓慢变化（~18.6 年）直接对应升交点退行对月球南北范围的调制。

输出: data/obs_envelope_beta_date.png
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FRAME = 'j2000' if '--j2000' in sys.argv else 'date'
DATA_FILE = f'data/envelope_beta_{FRAME}.txt'
OUT_FILE  = f'data/obs_envelope_beta_{FRAME}.png'
NODAL_YR  = 18.61   # 交点退行周期（年）


def load(path):
    days, beta, kind = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            days.append(float(p[0]))
            beta.append(float(p[1]))
            kind.append(int(p[2]))
    return np.array(days), np.array(beta), np.array(kind)


def main():
    if not os.path.exists(DATA_FILE):
        sys.exit(f'[Error] {DATA_FILE} not found. Run: make envelope')

    days, beta, kind = load(DATA_FILE)
    yr = days / 365.25

    pos = kind ==  1   # β 极大（月球偏北最远）
    neg = kind == -1   # β 极小（月球偏南最远）

    fig, ax = plt.subplots(figsize=(14, 5))

    ax.plot(yr[pos], beta[pos], color='#d62728', lw=0.5, alpha=0.8,
            label='β max (northernmost reach each orbit)')
    ax.plot(yr[neg], beta[neg], color='#1f77b4', lw=0.5, alpha=0.8,
            label='β min (southernmost reach each orbit)')
    ax.fill_between(yr[pos], beta[pos], 0, alpha=0.08, color='#d62728')
    ax.fill_between(yr[neg], beta[neg], 0, alpha=0.08, color='#1f77b4')
    ax.axhline(0, color='gray', lw=0.6, ls='--', label='Ecliptic (β = 0)')

    # 标记约 18.6 年节点退行周期
    for k, yr0 in enumerate(np.arange(0, yr[-1], NODAL_YR)):
        lbl = f'~{NODAL_YR} yr cycle' if k == 0 else None
        ax.axvline(yr0, color='orange', lw=0.8, alpha=0.6, ls=':', label=lbl)

    ax.set_xlabel('Years since 1900-01-01')
    ax.set_ylabel('Ecliptic latitude β (deg)')
    frame_label = 'SOFA date ecliptic' if FRAME == 'date' else 'fixed J2000 ecliptic'
    ax.set_title('Lunar β amplitude envelope — nodal regression modulates N/S reach\n'
                 f'(DE440, 1900–2100, {frame_label}, orange lines mark ~18.6 yr nodal period)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f'Saved: {OUT_FILE}')


if __name__ == '__main__':
    main()
