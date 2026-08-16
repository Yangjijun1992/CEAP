# CEAP — PMT Afterpulse 随机耦合本底评估模拟框架

PMT afterpulse 在 1 μs 时间窗口内随机耦合产生的实验背景噪声评估框架，
面向 RELICS 类液氙探测器（CEvNS 搜索）。

## 目标

1. 评估 PMT afterpulse 在 1 μs 窗口内随机耦合产生的背景噪声。
2. 评估大信号（宇宙线 μ 子 S1）过后 afterpulse 随机耦合形成的本底水平。

## 快速开始

```bash
conda activate py12
export PYTHONPATH=$PWD/src

# 运行演示（三种典型 PMT afterpulse 分布情形，逐 PMT PDF）
python -m ceap.cli --config config/afterpulse_scenarios.yaml

# 10 万次 μ子大规模模拟（假信号率/PE/宽度/Δt 直方图）
python -m ceap.tools.muon_s1_sim_100k --n-muons 100000 --rate 10

# 绘制 afterpulse 2D PDF（归一化前后对比）
python -m ceap.tools.plot_app_pdfs

# 运行单元测试
python -m pytest tests/
```

## 模拟框架

配置驱动的可插拔模块化设计（SIM 组），核心代码与具体数值解耦：

- **配置驱动**：探测器 / 阈值 / PDF / μ 子参数全部由 YAML 配置注入，默认配置与用户覆盖递归合并。
- **可插拔 PDF**：统一接口 `BaseAfterpulsePDF`，支持外部 `.npz` 文件（`FileAfterpulsePDF`）、
  逐 PMT PDF（`PerPMTPDF`）与占位实现（`PlaceholderAfterpulsePDF`）。
- **可复现**：随机种子统一管理，同一配置可复现整条时间线与结果。

数据流：

```
YAML 配置 ──► Config.loader（默认 + 用户覆盖递归合并）
                  │
     ┌────────────┼──────────────┐
     ▼            ▼              ▼
  Detector    AfterpulsePDF   MuonModel
 (探测器参数)  (AP 组外部数据)  (μ子 S1 生成)
     │            │              │
     └────────────┴──────┬───────┘
                         ▼
  1. TimelineSimulator：注入物理本底 + 暗噪声 → Hit 列表
  2. μ子主脉冲：到临时间 + S1 PE → MainPulse
  3. AfterpulseGenerator：每主脉冲按 PDF 独立为各 PMT 抽样 (时间, 电荷)
  4. WindowScanner：1 μs 窗扫描 → 超阈值事例 + 死区排除
                         ▼
      输出指标：本底率 / PE 谱 / 重叠 → JSON + NPY
```

| 模块 | 文件 | 职责 |
|------|------|------|
| 配置加载 | `src/ceap/config/loader.py` | YAML 递归合并与 Config 访问 |
| afterpulse PDF | `src/ceap/models/afterpulse_pdf.py` | 二维 PDF 接口与三种实现 |
| μ 子模型 | `src/ceap/models/muon.py` | 到临时间 / S1 PE 生成（M 组） |
| 探测器参数 | `src/ceap/models/detector.py` | 探测器参数容器（D 组） |
| 时间线 | `src/ceap/simulation/timeline.py` | 物理本底 + 暗噪声注入（SIM-01） |
| afterpulse 生成 | `src/ceap/simulation/afterpulse_gen.py` | 逐主脉冲按 PDF 抽样（SIM-02） |
| 窗口扫描 | `src/ceap/simulation/window_scanner.py` | 1 μs 滑动窗、阈值判定、死区排除（SIM-03/05） |
| 主流程 | `src/ceap/simulation/runner.py` | 流程编排与输出 |

## 模拟方法

**afterpulse 事件生成**

- 线性归一化模型（AP-02）：`p_ap_total = APP / <事件面积>`，每主 PE 以 `p_ap_total`
  概率独立产生 afterpulse（二项抽样）。
- 对主脉冲总 PE 数，从二维 PDF（时间差 × 电荷，dt 0–60 μs、npe 0.5–30.5 PE）独立为
  各 PMT 抽样 (dt, npe)；`afterpulse.pmt_pdfs` 可配置逐 PMT PDF（AP-05）。
- 二维 PDF 由 `tools/app_pdf_builder.py` 从 `pmt_analysis` 原始散点数据生成（`data/app_pdf_*.npz`）。

**时间线模拟（SIM-01）**

- 暗噪声按 `dark_rate_hz_per_pmt` 泊松注入；物理本底按 `background_rate_hz` 注入主脉冲。
- μ 子到临时间 = 通量 × 接收度泊松抽样（M-02）；
  S1 PE = 能损（dE/dx × 径迹长度）× 光产额 + 光子统计涨落（M-04）。

**窗口扫描与假信号判定（SIM-03/05）**

- 1 μs 滑动窗口（步长可配，默认 0.5 μs），累加窗内总 PE，判定是否落入信号阈值区间
  `[threshold_pe_min, threshold_pe_max]`。
- 死区排除（S-04）：S1 后 `s1_after_ns` 区间、μ 子 S2 区间剔除。
- 统计输出：假信号率（Hz）、PE 谱、时间宽度、距 μ 子 Δt。

**大规模模拟（向量化批处理）**

10 万次 μ 子模拟（`tools/muon_s1_sim_100k.py`）按 500 μ 子/batch 向量化：

1. 二项抽样每 PMT 的 afterpulse 事件数（`binomial(per_pmt_pe, p_ap_total)`）；
2. 向量化从二维 PDF 抽样 (dt, npe)；
3. 按 50 ns bin 累加电荷谱，滑动窗口求和判定假信号（[100 ns, 5 μs] 内、1 μs 窗 ∈ [120, 240] PE）。

## 探测器参数

探测器参数在 `config/settings.yaml` 的 `detector` 段集中配置（D 组，标注 [TBD] 的为后期真实值）：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `n_pmts` / `pmt_list` | `0` / `[]` | PMT 总数 / 列表（D-01）[TBD] |
| `single_pe.quantum_efficiency` | `0.0` | 量子效率（D-02）[TBD] |
| `single_pe.collection_efficiency` | `0.0` | 收集效率（D-02）[TBD] |
| `single_pe.gain` | `0.0` | 增益（D-02）[TBD] |
| `target.active_mass_kg` | `0.0` | 液氙有效质量（D-03）[TBD] |
| `target.light_yield_pe_keV` | `0.0` | 光产额（每 keV PE 数）（D-03）[TBD] |
| `target.quenching` | `{er: 0, nr: 0}` | 电子 / 核反冲淬灭因子（D-03）[TBD] |
| `daq.sample_rate_mhz` | `1.0` | 采样率（D-04）[TBD] |
| `daq.waveform_length_us` | `0.0` | 波形长度（决定 S1 死区时长）[TBD] |
| `daq.dead_time_us` | `0.0` | DAQ 死区时间（D-04）[TBD] |
| `background_rates` | `[]` | 背景事例率总表（D-05）[TBD] |
| `environment.depth_mwe` / `temperature_k` | `0.0` / `0.0` | 环境深度 / 温度（D-06）[TBD] |

**当前分析使用的参考场景**（10 万次 μ 子模拟）：

| 参数 | 数值 |
|------|------|
| PMT 数量 | 56 只（均匀接收） |
| μ 子 S1 总 PE | 1,000,000 PE（120 MeV 沉积） |
| 每 PMT 主脉冲 PE | ~17,857 PE |
| PMT APP | 5%（LV2358 PDF，app≈5.3%） |
| 信号窗 / 阈值 | 1 μs / 120–240 PE（CEvNS 阈值窗） |
| 死区建议 | S1 后 ≥ 5 μs |

## 项目结构

```
CEAP/
├── config/            # YAML 配置（探测器参数、阈值、PDF 映射）
├── data/              # afterpulse 二维 PDF 数据（data/app_pdf_*.npz）
├── docs/              # 需求文档、任务清单、框架设计、分析报告、进展跟踪
├── src/ceap/
│   ├── config/        # 配置加载
│   ├── models/        # afterpulse PDF、μ子模型、探测器参数
│   ├── simulation/    # 时间线、afterpulse 生成、窗口扫描、主流程
│   ├── tools/         # PDF 构建、绘图、μ子专项模拟
│   └── analysis/      # 验证与系统误差分析（预留）
├── tests/             # 单元测试
└── output/            # 运行输出（模拟结果、图，不入库）
```

## 文档

| 文档 | 说明 |
|------|------|
| [docs/PMT_Afterpulse_Background_Requirements.md](docs/PMT_Afterpulse_Background_Requirements.md) | 需求文档 |
| [docs/PMT_Afterpulse_Development_Task_List.md](docs/PMT_Afterpulse_Development_Task_List.md) | 开发任务清单 |
| [docs/PMT_Afterpulse_Simulation_Framework.md](docs/PMT_Afterpulse_Simulation_Framework.md) | 框架设计 |
| [docs/PMT_Afterpulse_Development_Progress.md](docs/PMT_Afterpulse_Development_Progress.md) | 开发进展跟踪 |
| [docs/PMT_Afterpulse_Muon_S1_Analysis_Report.md](docs/PMT_Afterpulse_Muon_S1_Analysis_Report.md) | μ子 S1 afterpulse 分析报告 |
| [docs/PMT_Afterpulse_Muon_S1_100k_Analysis_Report.md](docs/PMT_Afterpulse_Muon_S1_100k_Analysis_Report.md) | 10 万次 μ子大规模模拟报告 |
| [docs/app_area_vs_delta_time.md](docs/app_area_vs_delta_time.md) | APP 原始数据说明（pmt_analysis 导出） |

## 关键结果（10 万次 μ子，μ子率 10 Hz）

| 指标 | 数值 |
|------|------|
| 每 μ子假信号数 | 10.15 |
| 假信号率（@10 Hz μ子率） | 101.5 Hz |
| 假信号 PE 均值 | 181.2 PE |
| 假信号时间宽度均值 | 0.95 μs |
| 假信号距 μ子 Δt 均值 | 4.25 μs |

详见 [100k 分析报告](docs/PMT_Afterpulse_Muon_S1_100k_Analysis_Report.md)。

## 依赖

- Python ≥ 3.9
- numpy, PyYAML
- scipy（可选）, matplotlib（绘图）

## License

© Yangjijun 2026
