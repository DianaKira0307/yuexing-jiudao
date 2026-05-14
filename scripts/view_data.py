"""读取 lunar_state_1900_2100.bin 并打印统计信息
无第三方依赖，仅使用 Python 标准库 struct
"""

import struct
import math
import sys
import os

BIN_FILE = "data/lunar_state_1900_2100.bin"
META_FILE = "data/lunar_state_1900_2100.meta"

RECORD_BYTES = 64  # 8 个 double × 8 字节


def read_meta():
    """读取并打印元数据文件"""
    if not os.path.exists(META_FILE):
        print(f"[Warning] {META_FILE} 不存在")
        return
    print("=" * 60)
    print("  元数据")
    print("=" * 60)
    with open(META_FILE) as f:
        print(f.read(), end="")


def read_bin():
    """读取二进制文件并打印统计信息"""
    if not os.path.exists(BIN_FILE):
        print(f"[Error] {BIN_FILE} 不存在")
        sys.exit(1)

    file_size = os.path.getsize(BIN_FILE)
    n = file_size // RECORD_BYTES
    print(f"  文件大小: {file_size} bytes")
    print(f"  采样点数: {n}")
    print(f"  (预期: {n * RECORD_BYTES} / {file_size} — 一致)\n")

    with open(BIN_FILE, "rb") as f:
        raw = f.read()

    # 遍历所有记录，计算统计量
    r_min = float("inf")
    r_max = -float("inf")
    r_sum = 0.0
    v_min = float("inf")
    v_max = -float("inf")

    first = None
    last = None

    for i in range(n):
        off = i * RECORD_BYTES
        tdb1, tdb2, x, y, z, vx, vy, vz = struct.unpack("8d", raw[off:off + RECORD_BYTES])

        if i == 0:
            first = (tdb1, tdb2)
        if i == n - 1:
            last = (tdb1, tdb2)

        r = math.sqrt(x * x + y * y + z * z)
        v = math.sqrt(vx * vx + vy * vy + vz * vz)

        if r < r_min:
            r_min = r
        if r > r_max:
            r_max = r
        r_sum += r
        if v < v_min:
            v_min = v
        if v > v_max:
            v_max = v

    r_mean = r_sum / n

    tdb_start = first[0] + first[1]
    tdb_end = last[0] + last[1]

    print("=" * 60)
    print("  统计信息")
    print("=" * 60)
    print(f"  首采样点 TDB (JD): {tdb_start:.6f}")
    print(f"  末采样点 TDB (JD): {tdb_end:.6f}")
    print(f"  TDB 时间跨度:       {tdb_end - tdb_start:.1f} 天")
    print(f"  总采样数:           {n}")
    print(f"  时间步长:           {(tdb_end - tdb_start) / (n - 1):.6f} 天")
    print()
    print(f"  地心距离 (km):")
    print(f"    最小值:  {r_min:.3f}")
    print(f"    最大值:  {r_max:.3f}")
    print(f"    平均值:  {r_mean:.3f}")
    print(f"    振幅:    {r_max - r_min:.3f}")
    print()
    print(f"  地心速度 (km/s):")
    print(f"    最小值:  {v_min:.6f}")
    print(f"    最大值:  {v_max:.6f}")
    print()

    # 打印前 3 条记录作为样例
    print("=" * 60)
    print("  前 3 条记录")
    print("=" * 60)
    print(f"  {'tdb1':>16s} {'tdb2':>16s} {'x(km)':>14s} {'y(km)':>14s} {'z(km)':>14s} "
          f"{'vx(km/s)':>12s} {'vy(km/s)':>12s} {'vz(km/s)':>12s}")
    for i in range(min(3, n)):
        off = i * RECORD_BYTES
        vals = struct.unpack("8d", raw[off:off + RECORD_BYTES])
        print(f"  {vals[0]:16.6f} {vals[1]:16.6f} {vals[2]:14.3f} {vals[3]:14.3f} "
              f"{vals[4]:14.3f} {vals[5]:12.6f} {vals[6]:12.6f} {vals[7]:12.6f}")


def main():
    read_meta()
    print()
    read_bin()


if __name__ == "__main__":
    main()
