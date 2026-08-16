# PMT After Pulse 模拟框架设计文档

> 依据《PMT_Afterpulse_Development_Task_List.md》开发任务清单，搭建的
> Python 模拟软件框架（SIM 组）。采用**可插拔模块 + 外部 YAML 配置**结构：
> 探测器参数与 after pulse 二维 PDF 均为后期预留接口，当前用占位实现，
> 连接真实数据后无需改动核心逻辑。

---

## 1. 框架定位与设计原则

- **技术栈**：Python 3.9+，依赖 `numpy`、`PyYAML`（可选 `scipy`）。
- **配置驱动**：所有探测器/阈值/PDF/μ 子参数在 YAML 配置中预留，核心代码不与具体数值耦合。
- **可插拔**：after pulse PDF 通过统一接口（`BaseAfterpulsePDF`）注入，支持外部文件或参数化实现。
- **可复现**：随机种子统一管理，同一配置可复现整条时间线与结果。
- **可运行**：即使 PDF 与探测器参数未提供，也能用占位实现跑通全流程并输出指标。

---

## 2. 代码结构

```
CEAP/
├── pyproject.toml               # 包配置与依赖
├── conftest.py                  # 测试路径配置
├── config/
│   ├── settings.yaml            # 默认配置（探测器/PDF/μ子 占位）
│   ├── demo.yaml                # 暗噪声+本底 演示
│   └── demo_muon.yaml           # μ子 afterpulse 演示
├── src/ceap/
│   ├── __init__.py
│   ├── __main__.py              # python -m ceap 入口
│   ├── cli.py                   # 命令行入口
│   ├── config/
│   │   ├── loader.py            # YAML 加载 + 递归合并 + Config 访问
│   ├── models/
│   │   ├── afterpulse_pdf.py    # afterpulse 二维 PDF 接口与实现
│   │   ├── muon.py              # μ 子 S1 模型
│   │   └── detector.py          # 探测器参数容器 (D组预留)
│   ├── simulation/
│   │   ├── timeline.py          # SIM-01 时间线（物理本底+暗噪声）
│   │   ├── afterpulse_gen.py    # SIM-02 afterpulse 生成器
│   │   ├── window_scanner.py    # SIM-03/05 窗口扫描+指标
│   │   └── runner.py            # 主流程编排
│   ├── analysis/                # V 组分析验证（预留）
│   └── utils/                   # 工具
├── data/                        # 输入数据（PDF 文件等）
├── output/                      # 运行结果输出
└── tests/                       # 单元测试
```

---

## 3. 数据流 / 处理流程

```
                    ┌────────────────────────────────────────────┐
  YAML 配置  ──────►│ Config.loader (默认+用户覆盖递归合并)       │
                    └───────┬────────────────────────────────────┘
                            │
                            ▼
   ┌─────────────┐   ┌───────────────┐   ┌──────────────┐
   │Detector(D组)│   │AfterpulsePDF  │   │ MuonModel(M组)│
   │  参数占位    │   │  (AP组外部数据)│   │  S1 生成      │
   └─────────────┘   └──────┬────────┘   └──────┬───────┘
                            │                   │
                            ▼                   ▼
   ┌────────────────────────────────────────────────────┐
   │ 1. TimelineSimulator (SIM-01)                       │
   │    注入物理本底 + 暗噪声 → Hit 列表                  │
   ├────────────────────────────────────────────────────┤
   │ 2. μ子主脉冲 (M组): 到临时间 + S1 PE → MainPulse     │
   ├────────────────────────────────────────────────────┤
   │ 3. AfterpulseGenerator (SIM-02)                     │
   │    每主脉冲按 PDF 生成 afterpulse → 追加 Hit        │
   ├────────────────────────────────────────────────────┤
   │ 4. WindowScanner (SIM-03)                           │
   │    1μs 窗口扫描 → 超阈值事例 + 死区排除             │
   └──────────────────┬─────────────────────────────────┘
                      ▼
         输出指标 (S-05): 本底率/PE谱/重叠 → JSON+NPY
```

**关键步骤说明**

| 步骤 | 模块 | 说明 |
|------|------|------|
| 时间线生成 | `TimelineSimulator` | 按 `dark_rate_hz_per_pmt` 泊松注入暗噪声；按 `background_rate_hz` 注入物理本底主脉冲(S-03) |
| μ 子主脉冲 | `MuonModel` | 到临时间=通量×接收度泊松抽样；S1 PE=能损×光产额+涨落 |
| after pulse | `AfterpulseGenerator` | 对每个主脉冲的 PE 数，按 PDF 独立为各 PMT 抽样 (时间,电荷) |
| 窗口扫描 | `WindowScanner` | 滑动/步进 1μs 窗口，累加窗内总 PE，判定 [阈值区间] 与死区排除(S-04)，统计噪声事例率与 PE 谱 |

---

## 4. 配置说明（`config/settings.yaml`）

所有字段均有默认值；用户复制为 `settings.user.yaml` 或传入自定义 YAML 覆盖，
会与默认递归合并。**标注 `[TBD]` 的为后期提供项**。

| 配置段 | 关键字段 | 说明 | 状态 |
|--------|----------|------|------|
| `simulation` | `duration_us`, `seed`, `timeline.*` | 时长/随机种子/注入率(SIM-01) | 可用 |
| `simulation.timeline` | `dark_rate_hz_per_pmt`, `background_rate_hz` | 暗噪声/物理本底率 | **[TBD] D-05** |
| `detector` | `n_pmts`, `pmt_list`, `single_pe`, `target`, `daq`, `environment` | 探测器参数 | **[TBD] D 组全部** |
| `afterpulse` | `pdf_file`, `params`, `normalization.p_ap_*` | afterpulse 二维 PDF | **[TBD] AP 组** |
| `muon` | `flux_per_m2_hz`, `active_area_m2`, `dedx_mev_cm`, `track_length_mean_cm`, `light_yield_pe_keV` | μ 子 S1 | **[TBD] M 组** |
| `signal_window` | `length_us`, `step_us`, `threshold_pe_*`, `dead_zone` | 信号窗/阈值/死区(S 组) | 可用 |
| `output` | `save_hits`, `save_windows` | 输出开关 | 可用 |

---

## 5. Afterpulse 二维 PDF（AP 组，已接入真实数据）

统一接口：`BaseAfterpulsePDF`。提供三种实现：

1. **`FileAfterpulsePDF`** —— 从 `.npz` 文件加载二维直方图（由
   `src/ceap/tools/app_pdf_builder.py` 从 `pmt_analysis` 原始散点数据生成）。
   npz 含：
   - `h`: 二维直方图，shape `(len(dt_edges)-1, len(npe_edges)-1)`
   - `dt_edges`: 时间差轴（ns），默认 0–60 μs、100 ns 步长
   - `npe_edges`: 电荷轴（1 PE 步长，默认 0.5–30.5 PE）
   - `p_ap_total`: 每主 PE 全程总 after pulse 概率（AP-03）
   - `p_ap_1us`: 每主 PE 在 1 μs 窗内概率（AP-03）
   - `app`: 原始 APP = sum(ap_area_pe)/sum(main_area_pe)
   - `pmt_id`: 来源 PMT
   归一化模型（AP-02 线性）：`p_ap_total = APP / <事件面积>`，
   每主 PE 以 `p_ap_total` 概率独立产生 after pulse 事件。

2. **`PerPMTPDF`** —— 逐 PMT PDF 集合（AP-05），由 `afterpulse.pmt_pdfs`
   配置 `{pmt_index: 文件}` 构建，未配置的 PMT 用 `pdf_file` 兜底。

3. **`PlaceholderAfterpulsePDF`** —— 不产生 after pulse，用于流程调试。

### 典型 PMT 数据（data/app_pdf_*.npz）

| PMT | APP | p_ap_total | p_ap_1us | 分布情形 |
|-----|-----|-----------|----------|----------|
| LV2305 | 0.102 (~10%) | 0.0387 | 0.0246 | 中 after pulse |
| LV2358 | 0.053 (~5%)  | 0.0219 | 0.0153 | 低 after pulse |

> 注：LV2229（APP≈20%）因统计量太小（仅 7.5 万 afterpulse，对比百万级）暂不纳入分析参考。
> AP 组 PDF 已接入，抽样统计 APP 与原始值一致（相对偏差 < 2%）。
> 若 PDF 依赖主脉冲 PE 数（AP-02 非线性），重写 `BaseAfterpulsePDF.sample()` 即可。

---

## 6. 运行方式

```bash
# 激活环境 (py12)
conda activate py12
export PYTHONPATH=$PWD/src

# 默认配置运行（占位，仅跑通流程）
python -m ceap.cli --config config/settings.yaml

# 暗噪声+物理本底 演示
python -m ceap.cli --config config/demo.yaml

# μ子 afterpulse 演示（含外部 PDF）
python -m ceap.cli --config config/demo_muon.yaml

# 典型 PMT afterpulse 分布情形演示（逐 PMT PDF）
python -m ceap.cli --config config/afterpulse_scenarios.yaml

# 指定 run_id
python -m ceap.cli -c config/demo.yaml --run-id myrun
```

### 输出

- `output/<run_id>_summary.json`：汇总指标（主脉冲数、扫描窗口数、触发事例数、本底率 Hz、PE 谱、PDF 描述）。
- `output/<run_id>_windows.npy`：超阈值窗口明细 (时间, 总PE, 是否含物理S1)。可选。

---

## 7. 里程碑与后续接入项（对照任务清单）

| 阶段 | 内容 | 对应任务 | 框架现状 |
|------|------|----------|----------|
| M1 | 接入探测器参数 | T-D-01..06 | 占位，config 预留 |
| M1 | 接入 afterpulse PDF 数据 | T-AP-01..06 | ✅ 两个典型 PMT 已接入 (data/) |
| M2 | 固化 μ 子模型 | T-M-01..05 | 占位 S1 生成器 |
| M3 | 模拟软件可用 | T-SIM-01..05, T-S-01..03 | ✅ 核心可用 |
| M4 | 验证与系统误差 | T-V-01..04 | 分析模块预留 |
| M5 | 集成交付 | T-INT-01..04 | 主流程已集成 |

---

## 8. 后续工作（TBD 未接入项）

- **探测器参数**（D-02 单 PE 响应、D-05 背景率）：填 Config 后由 `Detector`/`TimelineSimulator` 读取。
- **更多 PMT PDF 数据**：用 `src/ceap/tools/app_pdf_builder.py` 从 `pmt_analysis`
  仓库批量生成其他 PMT 的 PDF 至 `data/`，配置 `afterpulse.pmt_pdfs` 映射。
- **μ 子模型细化**（M-01..05）：真实通量/角分布/能损参数。
- **验证模块**（V 组）：`analysis/` 下实现数据对比与灵敏度分析。
- **统计精度控制**（SIM-05）：按 `target_uncertainty_pct` 自动决定模拟时长。

---

*本文档随需求与代码变更同步更新。*
