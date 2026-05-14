"""月行九道 — 黄道交点漂移可视化

读取 data/node_crossings.txt，绘制升/降交点黄经随时间的变化。
退行率约 −19.3°/年，对应 ~18.6 年完成一圈，是"月行九道"最直接的观测特征。

输出: data/node_crossings.png
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_FILE = 'data/node_crossings.txt'
OUT_FILE  = 'data/obs_node_crossings.png'


def load(path):
    da, la, dd, ld = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            d, l, t = float(p[0]), float(p[1]), p[2]
            if t == 'A':
                da.append(d); la.append(l)
            else:
                dd.append(d); ld.append(l)
    return (np.array(da), np.array(la),
            np.array(dd), np.array(ld))


def unwrap_node(lam_deg):
    """对退行的交点黄经去卷绕（交点每步减少约 1.4°）"""
    rad = np.deg2rad(lam_deg)
    # 节点退行，相邻差期望为负值（~−0.025 rad/draconitic month）
    # unwrap 默认假设差值在 (−π, π]，对减量序列同样适用
    unwrapped = np.unwrap(rad)
    return np.rad2deg(unwrapped)


def main():
    if not os.path.exists(DATA_FILE):
        sys.exit(f'[Error] {DATA_FILE} not found. Run: make nodes')

    da, la, dd, ld = load(DATA_FILE)
    yr_a = da / 365.25
    yr_d = dd / 365.25

    la_uw = unwrap_node(la)
    ld_uw = unwrap_node(ld)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle('Lunar node drift — nodal regression (Nine Paths of the Moon)\n'
                 '(DE440, 1900–2100)', fontsize=12)

    # ---- 上行：散点（λ mod 360°），展示退行的"斜条纹"视觉 ----
    for ax, (yr, lam, label, color) in zip(axes[0], [
        (yr_a, la, 'Ascending node  (β: − → +)', '#d62728'),
        (yr_d, ld, 'Descending node (β: + → −)', '#1f77b4'),
    ]):
        ax.scatter(yr, lam, s=0.8, color=color, alpha=0.4, rasterized=True)
        ax.set_ylim(0, 360)
        ax.set_ylabel('Node longitude λ (deg, mod 360°)')
        ax.set_title(label + '  [mod 360°]')
        ax.grid(True, alpha=0.2)

    # ---- 下行：展开后的退行曲线 + 线性拟合 ----
    for ax, (yr, lam_uw, label, color) in zip(axes[1], [
        (yr_a, la_uw, 'Ascending node (unwrapped)', '#d62728'),
        (yr_d, ld_uw, 'Descending node (unwrapped)', '#1f77b4'),
    ]):
        ax.plot(yr, lam_uw, color=color, lw=0.5, alpha=0.7)

        p = np.polyfit(yr, lam_uw, 1)
        period = abs(360.0 / p[0])
        fit_line = np.polyval(p, yr[[0, -1]])
        ax.plot(yr[[0, -1]], fit_line, 'k--', lw=1.5,
                label=f'Regression: {p[0]:.2f} °/yr  →  period {period:.1f} yr')

        ax.set_xlabel('Years since 1900-01-01')
        ax.set_ylabel('Node longitude λ (deg, unwrapped)')
        ax.set_title(label)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f'Saved: {OUT_FILE}')


if __name__ == '__main__':
    main()
