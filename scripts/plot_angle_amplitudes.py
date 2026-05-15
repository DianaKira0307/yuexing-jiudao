"""月行九道 — 角变量实际振幅谱

从 data/lunar_state_1900_2100.bin 重建固定 J2000 黄道坐标下的
Δλ 和 β 序列，绘制 Hanning 相干增益校正后的单边半振幅谱。

输出: data/spectrum_angles_amplitude.png
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

BIN_FILE = 'data/lunar_state_1900_2100.bin'
OUT_FILE = 'data/spectrum_angles_amplitude.png'
N_FFT = 131072
P_MAX = 2.7 * 365.25
N_LABEL = 8
THRESHOLDS_DEG = (0.5, 0.1)
OBLIQUITY_J2000 = 0.409092804
RAD2DEG = 180.0 / np.pi


def load_angle_series():
    if not os.path.exists(BIN_FILE):
        sys.exit(f'[Error] {BIN_FILE} not found. Run: make run')

    raw = np.fromfile(BIN_FILE, dtype=np.float64)
    if raw.size % 8 != 0:
        sys.exit(f'[Error] {BIN_FILE} size is not a whole number of records')

    data = raw.reshape((-1, 8))
    x = data[:, 2]
    y = data[:, 3]
    z = data[:, 4]

    c = np.cos(OBLIQUITY_J2000)
    s = np.sin(OBLIQUITY_J2000)
    x_e = x
    y_e = c * y + s * z
    z_e = -s * y + c * z

    lam = np.arctan2(y_e, x_e)
    lam = np.where(lam < 0.0, lam + 2.0 * np.pi, lam)
    beta = np.arctan2(z_e, np.sqrt(x_e * x_e + y_e * y_e))

    t = np.arange(lam.size, dtype=np.float64)
    lam_unwrapped = np.unwrap(lam)
    slope, intercept = np.polyfit(t, lam_unwrapped, 1)
    dlambda = lam_unwrapped - (slope * t + intercept)

    return {
        'dlambda': dlambda * RAD2DEG,
        'beta': beta * RAD2DEG,
    }


def amplitude_spectrum(x):
    n = x.size
    if n > N_FFT:
        sys.exit(f'[Error] N_FFT={N_FFT} is shorter than the signal length {n}')

    x0 = x - np.mean(x)
    window = np.hanning(n)
    xw = x0 * window

    buf = np.zeros(N_FFT, dtype=np.float64)
    buf[:n] = xw

    spec = np.fft.rfft(buf)
    freq = np.fft.rfftfreq(N_FFT, d=1.0)

    amp = np.abs(spec) * 2.0 / np.sum(window)
    amp[0] = np.abs(spec[0]) / np.sum(window)
    if N_FFT % 2 == 0:
        amp[-1] = np.abs(spec[-1]) / np.sum(window)

    return freq, amp


def find_local_peaks(period, amp):
    peaks = []
    for i in range(1, amp.size - 1):
        if not (6.0 <= period[i] <= P_MAX):
            continue
        if amp[i] > amp[i - 1] and amp[i] > amp[i + 1]:
            peaks.append((amp[i], period[i], i))
    peaks.sort(reverse=True)
    return peaks


def fmt_period(days):
    if days < 200:
        return f'{days:.1f}d'
    if days < 3650:
        return f'{days / 365.25:.2f}yr'
    return f'{days / 365.25:.1f}yr'


def _xfmt(x, _):
    if x < 365:
        return f'{x:.0f}d'
    if x < 3650:
        return f'{x / 365.25:.1f}yr'
    return f'{x / 365.25:.0f}yr'


def draw_panel(ax, name, series, color, title):
    freq, amp = amplitude_spectrum(series)
    period = np.full_like(freq, np.inf)
    period[1:] = 1.0 / freq[1:]

    mask = (period >= 6.0) & (period <= P_MAX)
    ax.plot(period[mask], amp[mask], color=color, lw=0.7)

    peaks = find_local_peaks(period, amp)
    for a, p, i in peaks[:N_LABEL]:
        ax.scatter([p], [a], s=18, color='#cc2222', zorder=4)
        ax.annotate(
            f'{fmt_period(p)}\n{a:.2f}°',
            xy=(p, a),
            xytext=(5, 6),
            textcoords='offset points',
            fontsize=7.5,
            color='#aa1111',
            ha='left',
            va='bottom',
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(6, P_MAX)
    visible = amp[mask]
    positive = visible[visible > 0.0]
    if positive.size:
        ymin = max(1.0e-4, positive.max() * 1.0e-5)
        ymax = positive.max() * 2.5
        ax.set_ylim(ymin, ymax)
    for threshold in THRESHOLDS_DEG:
        ax.axhline(threshold, color='#444444', lw=1.0, ls='--', alpha=0.8)
        ax.text(
            np.sqrt(6.0 * P_MAX),
            threshold * 1.12,
            f'Threshold: {threshold:.1f} deg',
            ha='center',
            va='bottom',
            fontsize=8.5,
            color='#333333',
        )
    ax.set_xlabel('Period')
    ax.set_ylabel('Half-amplitude (deg)')
    ax.set_title(title)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_xfmt))
    ax.grid(True, which='major', ls='-', alpha=0.16)
    ax.grid(True, which='minor', ls='--', alpha=0.08)


def main():
    series = load_angle_series()

    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5))
    fig.subplots_adjust(hspace=0.35)
    fig.suptitle(
        f'Lunar angular amplitude spectra (DE440, fixed J2000 ecliptic, samples={len(series["beta"])}, FFT={N_FFT})',
        fontsize=11,
    )

    draw_panel(
        axes[0],
        'dlambda',
        series['dlambda'],
        '#2ca02c',
        'Longitude residual Δλ',
    )
    draw_panel(
        axes[1],
        'beta',
        series['beta'],
        '#d62728',
        'Ecliptic latitude β',
    )

    fig.savefig(OUT_FILE, dpi=150, bbox_inches='tight')
    print(f'Saved: {OUT_FILE}')


if __name__ == '__main__':
    main()
