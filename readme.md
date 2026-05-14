# 月行九道 (Lunar Nine Paths)

> 平台：macOS + gfortran  
> 目标：从 JPL DE440 星历提取月球地心状态矢量，通过 FFT 频谱分析研究月球轨道周期性，呼应古代"月行九道"概念。

---

## 一、依赖

| 依赖 | 说明 |
|------|------|
| gfortran ≥ 10 | Fortran 编译器 |
| IAU SOFA 静态库 | 已预编译为 `lib/libsofa.a`；如需自行编译见 [SOFA 官网](https://www.iausofa.org/) |
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

JPL DE440 星历文件（98 MB）不含于本仓库，需自行从 JPL 官方渠道获取 1-byte 记录版本，放置于项目根目录并命名为 `JPLEPH-DE440-1Bytes`。

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
├── data/          程序输出（.bin/.txt/.png 均在 .gitignore 中排除）
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
| `orbital_elements.f90` | 坐标变换 ICRF → 黄道球坐标 (r, λ, β) |
| `spectrum.f90` | 加窗、功率谱、寻峰 |

**`app/`——主程序**

| 文件 | 说明 |
|------|------|
| `main.f90` | 输出 `data/lunar_state_1900_2100.bin` |
| `analyze.f90` | 输出 `data/peaks_*.txt` |
| `envelope.f90` | 扩充：β 振幅包络 |
| `trajectory.f90` | 扩充：(λ, β) 轨迹 |
| `nodes.f90` | 扩充：黄道交点漂移 |

**`tests/`——单元测试**

| 文件 | 测试对象 |
|------|------|
| `test_fft.f90` | `mod_fft` |
| `test_lunar_extract.f90` | `mod_lunar_extract` |
| `test_orbital_elements.f90` | `mod_orbital_elements` |
| `test_spectrum.f90` | `mod_spectrum` |

| Python 脚本 | 说明 |
|------|------|
| `scripts/plot_spectra_mpl.py` | 绘制功率谱图 |
| `scripts/view_data.py` | 读取并统计状态矢量二进制文件 |
| `scripts/plot_envelope.py` | 绘制 β 振幅包络图 |
| `scripts/plot_trajectory.py` | 绘制黄道轨迹图 |
| `scripts/plot_nodes.py` | 绘制交点漂移图 |

---

## 三、模块依赖关系

```
mod_jpleph
mod_precision / mod_constants  ← mod_jpleph
mod_time_system                ← mod_precision  (链接 SOFA)
mod_lunar_extract              ← mod_precision, mod_jpleph, mod_time_system
mod_fft                        ← mod_precision
mod_orbital_elements           ← mod_precision
mod_spectrum                   ← mod_precision, mod_fft
```

---

## 四、构建与运行

```bash
make              # 查看帮助
make build        # 编译 extract_main
make run          # 生成 data/lunar_state_1900_2100.bin
make analyze      # 生成 data/peaks_*.txt
make plot         # 绘制功率谱图
make view         # 打印状态矢量统计

# 扩充分析（九道可视化）
make observe      # 生成 envelope/trajectory/nodes 数据
make observe_plot # 生成三张图
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

## 五、主要频谱结果

| 变量 | 主峰周期 | 物理含义 |
|------|---------|---------|
| β（黄纬） | 27.2 d | 交点月 |
| β | 5.98 yr | 升交点退行相关周期 |
| Δλ（黄经速率） | 27.6 d | 近点月 |
| Δλ | 366.1 d | 年差项 |
| r（地心距） | 27.6 d | 近点月 |
| r | 206.1 d | ≈ 7 朔望月（食年相关） |

---

## 六、声明

**AI 辅助说明**

本项目的代码编写与调试过程中使用了 [Claude](https://claude.ai)（Anthropic）与 [DeepSeek](https://www.deepseek.com) 作为辅助工具，用于代码生成、错误排查与方案讨论。目前呈现的代码均由 AI 生成，暂未进行复核。项目的封面图片由 ChatGPT 生成。

**代码来源说明**

`src/jpleph.f90`（JPL 星历文件读取）与 `src/time_and_coord.f90`（时间系统转换，封装 IAU SOFA）中的核心逻辑基于已有的参考实现改写而来，并根据本项目需求进行了适配。

---

*月行九道项目 | 2026 年 5 月*
