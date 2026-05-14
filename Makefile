# Makefile for Lunar Nine Paths Phase 1
# 月行九道 Phase 1: 月球地心状态矢量提取
#
# Targets:
#   make build  — 编译 extract_main
#   make run    — 编译并运行 extract_main (生成 data/*.bin + data/*.meta)
#   make test   — 编译并运行 test_lunar_extract (生成 data/test_lunar_extract.log)
#   make clean  — 清理 build 目录

FC       = gfortran
FFLAGS   = -O2 -std=f2018 -Wall -fimplicit-none -J build
SOFA_LIB ?= lib/libsofa.a
BUILD_DIR = build
DATA_DIR  = data
PYTHON    = .venv/bin/python3
PY_DIR    = scripts
VENV_SENTINEL = .venv/.installed

# 本项目使用的源文件（按编译顺序）
MOD_SRCS = jpleph.f90 utilities.f90 time_and_coord.f90 lunar_extract.f90 fft.f90 orbital_elements.f90 spectrum.f90
MOD_OBJS = $(addprefix $(BUILD_DIR)/, $(MOD_SRCS:.f90=.o))

.PHONY: all build run test test_fft test_orbital_elements test_spectrum analyze plot plot_svg view \
        envelope trajectory nodes plot_envelope plot_trajectory plot_nodes \
        observe observe_plot deploy venv clean clean-venv help

help:
	@echo "Targets: build | run | test | test_fft | analyze | plot | view | venv | clean"
	@echo "  build    : 编译 extract_main"
	@echo "  run      : 编译并运行 extract_main (生成 data/*.bin + data/*.meta)"
	@echo "  test     : 编译并运行 test_lunar_extract (生成 data/test_lunar_extract.log)"
	@echo "  test_fft : 编译并运行 test_fft (生成 data/test_fft.log)"
	@echo "  test_orbital_elements : 编译并运行 test_orbital_elements (生成 data/test_orbital_elements.log)"
	@echo "  test_spectrum : 编译并运行 test_spectrum (生成 data/test_spectrum.log)"
	@echo "  analyze  : 编译并运行 lunar_analyze (生成 data/peaks_*.txt)"
	@echo "  plot     : 运行 plot_spectra.py (生成 data/*.svg)"
	@echo "  view     : 运行 view_data.py (打印统计信息)"
	@echo "--- 扩充分析 ---"
	@echo "  envelope      : 编译并运行 lunar_envelope (生成 data/envelope_beta.txt)"
	@echo "  trajectory    : 编译并运行 lunar_trajectory (生成 data/traj_ecliptic.txt)"
	@echo "  nodes         : 编译并运行 lunar_nodes (生成 data/node_crossings.txt)"
	@echo "  plot_envelope : 绘制 β 振幅包络图 (data/envelope_beta.png)"
	@echo "  plot_trajectory: 绘制黄道轨迹图 (data/traj_ecliptic.png)"
	@echo "  plot_nodes    : 绘制交点漂移图 (data/node_crossings.png)"
	@echo "  observe       : 依次运行三个计算程序"
	@echo "  observe_plot  : 依次生成三张图"
	@echo "  deploy        : 一键全流程 (run+analyze+observe+所有绘图)"
	@echo "  venv     : 创建/更新 Python 虚拟环境（.venv/）"
	@echo "  clean    : 清理 build 目录"
	@echo "  clean-venv : 删除 .venv/"

# --- Python 虚拟环境 ---
$(VENV_SENTINEL): requirements.txt
	python3 -m venv .venv
	.venv/bin/pip install --quiet -r requirements.txt
	@touch $(VENV_SENTINEL)

venv: $(VENV_SENTINEL)

# --- 可执行文件 ---
build: $(BUILD_DIR)/extract_main

$(BUILD_DIR)/extract_main: $(MOD_OBJS) $(BUILD_DIR)/main.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/test_lunar_extract: $(MOD_OBJS) $(BUILD_DIR)/test_lunar_extract.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/test_fft: $(MOD_OBJS) $(BUILD_DIR)/test_fft.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/test_orbital_elements: $(MOD_OBJS) $(BUILD_DIR)/test_orbital_elements.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/test_spectrum: $(MOD_OBJS) $(BUILD_DIR)/test_spectrum.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/lunar_analyze: $(MOD_OBJS) $(BUILD_DIR)/analyze.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

# --- 通用编译规则 ---
$(BUILD_DIR)/%.o: src/%.f90
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: app/%.f90
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: tests/%.f90
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(DATA_DIR):
	mkdir -p $(DATA_DIR)

# --- 运行 ---
run: $(BUILD_DIR)/extract_main | $(DATA_DIR)
	./$(BUILD_DIR)/extract_main

test: $(BUILD_DIR)/test_lunar_extract | $(DATA_DIR)
	./$(BUILD_DIR)/test_lunar_extract

test_fft: $(BUILD_DIR)/test_fft | $(DATA_DIR)
	./$(BUILD_DIR)/test_fft

test_orbital_elements: $(BUILD_DIR)/test_orbital_elements | $(DATA_DIR)
	./$(BUILD_DIR)/test_orbital_elements

test_spectrum: $(BUILD_DIR)/test_spectrum | $(DATA_DIR)
	./$(BUILD_DIR)/test_spectrum

analyze: $(BUILD_DIR)/lunar_analyze | $(DATA_DIR)
	./$(BUILD_DIR)/lunar_analyze

plot: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/plot_spectra_mpl.py

plot_svg: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/plot_spectra.py --svg

view: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/view_data.py

# --- 扩充分析：可执行文件 ---
$(BUILD_DIR)/lunar_envelope: $(MOD_OBJS) $(BUILD_DIR)/envelope.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/lunar_trajectory: $(MOD_OBJS) $(BUILD_DIR)/trajectory.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

$(BUILD_DIR)/lunar_nodes: $(MOD_OBJS) $(BUILD_DIR)/nodes.o
	$(FC) $(FFLAGS) $^ "$(SOFA_LIB)" -o $@

# --- 扩充分析：计算 ---
envelope: $(BUILD_DIR)/lunar_envelope | $(DATA_DIR)
	./$(BUILD_DIR)/lunar_envelope

trajectory: $(BUILD_DIR)/lunar_trajectory | $(DATA_DIR)
	./$(BUILD_DIR)/lunar_trajectory

nodes: $(BUILD_DIR)/lunar_nodes | $(DATA_DIR)
	./$(BUILD_DIR)/lunar_nodes

# --- 扩充分析：绘图 ---
plot_envelope: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/plot_envelope.py

plot_trajectory: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/plot_trajectory.py

plot_nodes: $(VENV_SENTINEL) | $(DATA_DIR)
	$(PYTHON) $(PY_DIR)/plot_nodes.py

# --- 扩充分析：组合目标 ---
observe: envelope trajectory nodes

observe_plot: plot_envelope plot_trajectory plot_nodes

# --- 一键全流程部署 ---
deploy: run analyze observe plot observe_plot
	@echo "=== Deploy complete. All figures saved to data/ ==="

# --- 模块间依赖（确保 .mod 文件按正确顺序生成）---
# jpleph 无模块依赖 → 最先编译
# utilities(mod_constants)  use mod_jpleph → jpleph.o 必须先编译
# time_and_coord            use mod_precision → utilities.o 必须先编译
# lunar_extract             use mod_precision, mod_jpleph, mod_time_system
$(BUILD_DIR)/utilities.o: $(BUILD_DIR)/jpleph.o
$(BUILD_DIR)/time_and_coord.o: $(BUILD_DIR)/utilities.o
$(BUILD_DIR)/lunar_extract.o: $(BUILD_DIR)/utilities.o $(BUILD_DIR)/jpleph.o $(BUILD_DIR)/time_and_coord.o
$(BUILD_DIR)/fft.o: $(BUILD_DIR)/utilities.o
$(BUILD_DIR)/main.o: $(BUILD_DIR)/lunar_extract.o
$(BUILD_DIR)/test_lunar_extract.o: $(BUILD_DIR)/lunar_extract.o
$(BUILD_DIR)/orbital_elements.o: $(BUILD_DIR)/utilities.o
$(BUILD_DIR)/spectrum.o: $(BUILD_DIR)/fft.o
$(BUILD_DIR)/test_fft.o: $(BUILD_DIR)/fft.o
$(BUILD_DIR)/test_orbital_elements.o: $(BUILD_DIR)/orbital_elements.o $(BUILD_DIR)/lunar_extract.o
$(BUILD_DIR)/test_spectrum.o: $(BUILD_DIR)/spectrum.o
$(BUILD_DIR)/analyze.o:     $(BUILD_DIR)/spectrum.o $(BUILD_DIR)/lunar_extract.o $(BUILD_DIR)/orbital_elements.o
$(BUILD_DIR)/envelope.o:   $(BUILD_DIR)/orbital_elements.o $(BUILD_DIR)/utilities.o
$(BUILD_DIR)/trajectory.o: $(BUILD_DIR)/orbital_elements.o $(BUILD_DIR)/utilities.o
$(BUILD_DIR)/nodes.o:      $(BUILD_DIR)/orbital_elements.o $(BUILD_DIR)/utilities.o

# --- SOFA 库缺失时的提示 ---
$(SOFA_LIB):
	@echo "Missing SOFA library: $(SOFA_LIB)"
	@echo "Build it with: make -C ../../../lib/iau-sofa/src"
	@exit 1

# --- 清理 ---
clean:
	rm -rf $(BUILD_DIR)

clean-venv:
	rm -rf .venv
