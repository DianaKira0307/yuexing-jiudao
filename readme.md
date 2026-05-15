# 月行九道 (Lunar Nine Paths)

> 平台：macOS + gfortran  
> 目标：从 JPL DE440 星历提取月球地心状态矢量，通过频谱分析、黄纬包络、黄道轨迹和交点退行验证“月行九道”的几何与周期特征。

---

## 一、依赖

| 依赖 | 说明 |
|------|------|
| gfortran >= 10 | Fortran 编译器 |
| IAU SOFA 静态库 | 已预编译为 `lib/libsofa.a`；时间转换和 SOFA date-ecliptic 坐标转换会链接它 |
| JPL DE440 星历 | **不含于仓库**，需自行从 JPL 官方渠道获取，见下文说明 |
| Python 3 + NumPy + Matplotlib | 用于绘图脚本（`scripts/`），依赖托管在项目根目录 `.venv/` |

### Python 环境初始化

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

之后 `make plot` / `make observe_plot` 会自动调用 `.venv/bin/python3`。

### 获取 DE440 星历

JPL DE440 星历文件（约 98 MB）不含于本仓库，需自行从 JPL 官方渠道获取 1-byte 记录版本，放置于项目根目录并命名为：

```text
JPLEPH-DE440-1Bytes
```

---

## 二、文件结构

```
moon-orbit-analyze/
├── src/           Fortran 模块（库代码）
├── app/           Fortran 主程序（可执行入口）
├── tests/         Fortran 单元测试
├── scripts/       Python 可视化脚本
├── lib/
│   └── libsofa.a  IAU SOFA 静态库
├── data/          程序输出（.txt/.png 入仓；.bin/.meta 忽略）
├── build/         编译中间文件（gitignore）
├── Makefile
└── readme.md
```

**`src/`——模块**

| 文件 | 说明 |
|------|------|
| `jpleph.f90` | JPL 星历读取库封装 (`mod_jpleph`) |
| `utilities.f90` | 精度定义 + 物理常数 |
| `time_and_coord.f90` | 时间系统转换（封装 IAU SOFA） |
| `lunar_extract.f90` | 批量提取月球地心状态矢量 |
| `fft.f90` | Radix-2 Cooley-Tukey 复数 FFT |
| `orbital_elements.f90` | 坐标变换、黄道球坐标、轨道根数、Poincare 变量 |
| `spectrum.f90` | 加窗、功率谱、寻峰 |

**`app/`——主程序**

| 文件 | 说明 |
|------|------|
| `main.f90` | 输出 `data/lunar_state_1900_2100.bin/.meta` |
| `analyze.f90` | 频谱分析，输出 `data/peaks_{r,dlambda,beta}.txt` |
| `envelope.f90` | β 振幅包络，输出 J2000/date 两套结果 |
| `trajectory.f90` | `(λ, β)` 轨迹，输出 J2000/date 两套结果 |
| `nodes.f90` | 黄道交点漂移，输出 J2000/date 两套结果 |

**`tests/`——单元测试**

| 文件 | 测试对象 |
|------|------|
| `test_fft.f90` | `mod_fft` |
| `test_lunar_extract.f90` | `mod_lunar_extract` |
| `test_orbital_elements.f90` | `mod_orbital_elements` |
| `test_spectrum.f90` | `mod_spectrum` |

**`scripts/`——可视化**

| Python 脚本 | 说明 |
|------|------|
| `scripts/plot_spectra_mpl.py` | 绘制 `r / Δλ / β` 功率谱图 |
| `scripts/plot_angle_amplitudes.py` | 绘制 `Δλ / β` 实际角度半振幅谱 |
| `scripts/view_data.py` | 读取并统计状态矢量二进制文件 |
| `scripts/plot_envelope.py` | 绘制 β 振幅包络图，默认读取 date-ecliptic，`--j2000` 切换固定 J2000 |
| `scripts/plot_trajectory.py` | 绘制黄道轨迹图，默认读取 date-ecliptic，`--j2000` 切换固定 J2000 |
| `scripts/plot_nodes.py` | 绘制交点漂移图，默认读取 date-ecliptic，`--j2000` 切换固定 J2000 |

---

## 三、模块依赖关系

```
mod_jpleph
mod_precision / mod_constants  <- mod_jpleph
mod_time_system                <- mod_precision  (链接 SOFA 时间系统)
mod_lunar_extract              <- mod_precision, mod_jpleph, mod_time_system
mod_fft                        <- mod_precision
mod_orbital_elements           <- mod_precision  (SOFA iau_ECM06 用于 date-ecliptic)
mod_spectrum                   <- mod_precision, mod_fft
```

---

## 四、坐标系约定

项目现在区分两类黄道坐标系：

| 模式 | 文件后缀 | 用途 |
|------|----------|------|
| 固定 J2000 黄道近似 | `_j2000` | 动力学/统一参考系比较；使用固定 J2000 黄赤交角旋转 |
| SOFA 历元黄道坐标系 | `_date` | 模拟当时观测语义；使用 `iau_ECM06(date1,date2,RM)` 转到 mean ecliptic/equinox of date |

频谱分析 `make analyze` 保持使用固定 J2000 黄道近似，避免把坐标系随时间变化混入动力学周期。

观测类分析 `make observe` 同时输出两套文本结果：

```text
data/envelope_beta_j2000.txt
data/envelope_beta_date.txt
data/traj_ecliptic_j2000.txt
data/traj_ecliptic_date.txt
data/node_crossings_j2000.txt
data/node_crossings_date.txt
```

绘图脚本默认使用 `_date` 文件，生成：

```text
data/obs_envelope_beta_date.png
data/obs_traj_ecliptic_date.png
data/obs_node_crossings_date.png
```

如需固定 J2000 图：

```bash
.venv/bin/python3 scripts/plot_envelope.py --j2000
.venv/bin/python3 scripts/plot_trajectory.py --j2000
.venv/bin/python3 scripts/plot_nodes.py --j2000
```

---

## 五、构建与运行

```bash
make              # 查看帮助
make build        # 编译 extract_main
make run          # 生成 data/lunar_state_1900_2100.bin/.meta
make analyze      # 生成 data/peaks_{r,dlambda,beta}.txt
make plot         # 绘制频谱图 data/spectrum_ecliptic.png
make plot_angles  # 绘制 Δλ/β 实际角度半振幅谱
make view         # 打印状态矢量统计

# 扩充分析（月行九道观测可视化）
make observe      # 生成 envelope/trajectory/nodes 的 _j2000 与 _date 文本数据
make observe_plot # 默认绘制 date-ecliptic 三张图
make deploy       # 一键全流程

# 单元测试
make test
make test_fft
make test_orbital_elements
make test_spectrum

make clean        # 清理 build/
```

如需指定 SOFA 库路径：

```bash
make SOFA_LIB=/path/to/libsofa.a build
```

---

## 六、主要输出文件

| 文件 | 说明 | 是否入仓 |
|------|------|----------|
| `data/lunar_state_1900_2100.bin` | DE440 月球地心状态矢量二进制文件 | 否 |
| `data/lunar_state_1900_2100.meta` | 二进制文件元数据 | 否 |
| `data/peaks_r.txt` | 地心距 `r` 的谱峰表 | 是 |
| `data/peaks_dlambda.txt` | 黄经残差 `Δλ` 的谱峰表，含 `amp_deg` 半振幅列 | 是 |
| `data/peaks_beta.txt` | 黄纬 `β` 的谱峰表，含 `amp_deg` 半振幅列 | 是 |
| `data/envelope_beta_{j2000,date}.txt` | 黄纬包络极值表 | 是 |
| `data/traj_ecliptic_{j2000,date}.txt` | 黄道轨迹 `(λ, β)` 表 | 是 |
| `data/node_crossings_{j2000,date}.txt` | 黄道交点穿越表 | 是 |
| `data/spectrum_angles_amplitude.png` | `Δλ / β` 实际角度半振幅谱 | 是 |
| `data/*.png` | 频谱图、包络图、轨迹图、交点漂移图 | 是 |

---

## 七、主要频谱结果

**采样参数**：DE440 实际采样点数 `73050`，采样率 1 次/天，覆盖 1900-01-01 至 2100-01-01。

**FFT 参数**：实际样本不变，补零到 `131072` 点（2^17）执行 radix-2 FFT。谱峰表中：

```text
Samples: 73050
FFT length: 131072
```

| 变量 | 主峰周期 | 物理含义 |
|------|---------|---------|
| `β`（黄纬） | 27.21 d | 交点月相关的黄纬振荡 |
| `β` | 2184.53 d | 长周期调制项 |
| `Δλ`（黄经残差） | 27.55 d | 月球轨道主周期项 |
| `Δλ` | 365.10 d | 年差项 |
| `r`（地心距） | 27.55 d | 近点月相关的地心距振荡 |
| `r` | 31.81 d / 14.77 d | 月球轨道摄动和谐波项 |

`make plot` 输出的是相对功率谱，适合比较周期强弱；`make plot_angles` 输出的是角变量半振幅谱，纵轴单位为度，并标出 `0.5 deg` 与 `0.1 deg` 两条参考阈值线。当前主要角度峰约为：

| 变量 | 主峰周期 | 半振幅 |
|------|----------|--------|
| `Δλ` | 27.55 d | 约 6.18 deg |
| `β` | 27.21 d | 约 4.83 deg |

---

## 八、月行九道验证结果

项目通过四条证据链验证“月行九道”：

1. **频谱分析**：`r / Δλ / β` 的谱峰显示月球运动存在清晰的月周期、半月周期、年周期和长周期调制。
2. **黄纬包络**：`β` 的南北极值显示月球相对黄道的活动范围约为 ±5.3°。
3. **黄道轨迹**：`(λ, β)` 轨迹在 1 个月、1 年、9.3 年、18.6 年窗口中逐渐展开为黄道附近的路径带。
4. **交点退行**：`β=0` 的升/降交点持续退行，是“月行九道”最直接的几何证据。

当前结果：

| 坐标模式 | 升/降交点退行率 | 周期 |
|----------|----------------|------|
| 固定 J2000 黄道近似 | -19.3552 deg/yr | 18.5996 yr |
| SOFA 历元黄道坐标系 | -19.3413 deg/yr | 18.6130 yr |

两者差异来自坐标定义不同：固定 J2000 用统一参考平面，date-ecliptic 使用每个采样时刻的当日平黄道/春分点。

---

## 九、声明

**AI 辅助说明**

本项目的代码编写与调试过程中使用了 [Claude](https://claude.ai)（Anthropic）、[DeepSeek](https://www.deepseek.com) 与 ChatGPT/Codex 作为辅助工具，用于代码生成、错误排查、文档整理与方案讨论。项目的封面图片由 ChatGPT 生成。

**代码来源说明**

`src/jpleph.f90`（JPL 星历文件读取）与 `src/time_and_coord.f90`（时间系统转换，封装 IAU SOFA）中的核心逻辑基于已有的参考实现改写而来，并根据本项目需求进行了适配。`src/orbital_elements.f90` 中的 date-ecliptic 坐标转换使用 IAU SOFA `iau_ECM06`。

---

*月行九道项目 | 2026 年 5 月*
