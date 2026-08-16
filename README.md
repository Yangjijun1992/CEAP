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
