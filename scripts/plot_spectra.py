"""Lunar Nine Paths — 频谱可视化
零依赖：纯 Python 标准库，用 SVG 矢量图代替 matplotlib。

用法:
  python3 plot_spectra.py              # 仅打印文本统计
  python3 plot_spectra.py --svg        # 生成 SVG 图 (data/*.svg)
"""

import struct
import math
import sys
import os

BIN_FILE = "data/lunar_state_1900_2100.bin"
RECORD_BYTES = 64
VAR_NAMES = ["x", "y", "z", "vx", "vy", "vz"]
VAR_TITLES = [
    "X position (km)", "Y position (km)", "Z position (km)",
    "Vx velocity (km/s)", "Vy velocity (km/s)", "Vz velocity (km/s)",
]
N_DAYS_PLOT = 365 * 5


# ── 读二进制 ──

def read_binary(n_limit=65536):
    if not os.path.exists(BIN_FILE):
        print(f"[Error] {BIN_FILE} 不存在。先运行 make run")
        sys.exit(1)
    file_size = os.path.getsize(BIN_FILE)
    n_total = file_size // RECORD_BYTES
    n = min(n_total, n_limit)
    print(f"采样: {n}/{n_total}  跨度: {n / 365.25:.1f} 年\n")
    with open(BIN_FILE, "rb") as f:
        raw = f.read(n * RECORD_BYTES)
    tdb, data = [], [[] for _ in range(6)]
    for i in range(n):
        off = i * RECORD_BYTES
        vals = struct.unpack("8d", raw[off:off + 64])
        tdb.append(vals[0] + vals[1])
        for j in range(6):
            data[j].append(vals[2 + j])
    return tdb, data


# ── SVG 辅助 ──

def _svg_header(w, h, title):
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}"' \
           f' style="font-family:monospace;background:#fff">\n' \
           f'<text x="{w//2}" y="24" text-anchor="middle" font-size="16" ' \
           f'font-weight="bold">{_esc(title)}</text>\n'


def _esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _axes(x0, y0, w, h, xlabel, ylabel):
    """返回绘制坐标轴的 SVG 片段"""
    parts = [
        # 边框
        f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#333" stroke-width="1"/>',
        # 标签
        f'<text x="{x0 + w // 2}" y="{y0 + h + 36}" text-anchor="middle" font-size="12">{_esc(xlabel)}</text>',
        f'<text x="{x0 - 8}" y="{y0 + h // 2}" text-anchor="end" font-size="12" '
        f'transform="rotate(-90,{x0 - 8},{y0 + h // 2})">{_esc(ylabel)}</text>',
    ]
    return "\n".join(parts)


def _polyline(x0, y0, w, h, xs, ys, color="#3465a4", stroke_w=0.7):
    """绘制折线。xs, ys 是数据坐标列表，自动缩放到绘图区"""
    if not xs or not ys:
        return ""
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1
    pad = 0.05
    xr = xmax - xmin
    yr = ymax - ymin
    pts = []
    for xi, yi in zip(xs, ys):
        px = x0 + (xi - xmin) / xr * (w - 2 * pad * w) + pad * w
        py = y0 + h - (yi - ymin) / yr * (h - 2 * pad * h) - pad * h
        pts.append(f"{px:.1f},{py:.1f}")
    return f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="{stroke_w}"/>'


def _vline(x0, y0, h, x_frac, color="#888", style="dashed"):
    return f'<line x1="{x0 + x_frac:.1f}" y1="{y0}" x2="{x0 + x_frac:.1f}" ' \
           f'y2="{y0 + h}" stroke="{color}" stroke-dasharray="4,3" stroke-width="0.5"/>'


def _label(x0, y0, w, h, x_frac, y_frac, text, color="#333"):
    return f'<text x="{x0 + x_frac:.1f}" y="{y0 + h * (1 - y_frac) - 4:.1f}" ' \
           f'fill="{color}" font-size="10">{_esc(text)}</text>'


# ── 绘图函数 ──

def draw_timeseries(tdb, data, name, title, outpath):
    """时间序列折线图 → SVG"""
    n = min(N_DAYS_PLOT, len(tdb))
    days = [tdb[i] - tdb[0] for i in range(n)]
    W, H = 900, 240
    X0, Y0 = 70, 50
    PW, PH = 780, 160

    svg = _svg_header(W, 3 * H + 40, f"Lunar Geocentric {title} (first 5 years)")

    for row in range(3):
        y_off = Y0 + row * H
        base_idx = 0 if name == "position" else 3
        svg += _axes(X0, y_off, PW, PH, "Days since 1900-01-01", VAR_TITLES[base_idx + row])
        svg += _polyline(X0, y_off, PW, PH, days, data[row][:n])

    svg += "</svg>"
    with open(outpath, "w") as f:
        f.write(svg)
    print(f"  Saved: {outpath}")


def draw_spectrum(name, outpath):
    """谱峰图 → SVG"""
    path = f"data/peaks_{name}.txt"
    if not os.path.exists(path):
        print(f"  Skip: {path} not found")
        return

    freqs, powers, periods = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 5:
                freqs.append(float(parts[2]))
                periods.append(float(parts[3]))
                powers.append(float(parts[4]))

    if not freqs:
        return

    W, H = 900, 350
    X0, Y0 = 80, 50
    PW, PH = 760, 240

    # 对数功率
    log_pow = [math.log10(max(p, 1e-30)) for p in powers]
    max_n, min_n = max(powers), min(powers)
    ymin, ymax = math.log10(max(min_n, 1e-30)), math.log10(max_n)

    svg = _svg_header(W, H + 40, f"{name.upper()}: Power Spectrum (log scale)")

    # 框 + 标签
    svg += _axes(X0, Y0, PW, PH, "Frequency (cycles/day)", "log10(Power)")

    # 网格线（虚线）
    for frac in [0.2, 0.4, 0.6, 0.8]:
        svg += _vline(X0, Y0, PH, frac * PW, "#ddd")
        yl = Y0 + PH * (1 - frac)
        svg += f'<line x1="{X0}" y1="{yl}" x2="{X0 + PW}" y2="{yl}" stroke="#ddd" stroke-dasharray="2,2" stroke-width="0.5"/>'

    # 折线
    xs = freqs
    ys = log_pow
    svg += _polyline(X0, Y0, PW, PH, xs, ys, color="#3465a4", stroke_w=0.8)

    # 标注前 8 个峰值
    for j in range(min(8, len(freqs))):
        xf = (freqs[j] - min(xs)) / (max(xs) - min(xs) + 1e-30)
        yf = (log_pow[j] - ymin) / (ymax - ymin + 1e-30)
        px = X0 + xf * PW
        py = Y0 + PH * (1 - yf)
        label = f"{periods[j]:.1f}d" if periods[j] < 1000 else f"{periods[j] / 365.25:.2f}yr"
        svg += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3" fill="#cc0000"/>'
        svg += f'<text x="{px + 6:.1f}" y="{py - 4:.1f}" fill="#cc0000" font-size="9">{_esc(label)}</text>'

    svg += "</svg>"
    with open(outpath, "w") as f:
        f.write(svg)
    print(f"  Saved: {outpath}")


def print_peak_tables():
    print("=" * 60)
    print("  谱峰表")
    print("=" * 60)
    for name in VAR_NAMES:
        path = f"data/peaks_{name}.txt"
        if not os.path.exists(path):
            continue
        print(f"\n  --- {name} ---")
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) >= 5:
                        rank, idx, freq, period, power = parts[:5]
                        print(f"    #{rank:>3s}: f={freq} cyc/day  P={period:>8s} days  power={power}")
    print()


# ── 主入口 ──

def main():
    do_svg = "--svg" in sys.argv
    print("Lunar Nine Paths — Spectrum Analysis")
    print("=" * 60)
    tdb, data = read_binary()

    if do_svg:
        draw_timeseries(tdb, data[:3], "position", "Position", "data/timeseries_position.svg")
        draw_timeseries(tdb, data[3:], "velocity", "Velocity", "data/timeseries_velocity.svg")
        for name in ["x", "z"]:
            draw_spectrum(name, f"data/spectrum_{name}.svg")
        print()

    print_peak_tables()


if __name__ == "__main__":
    main()
