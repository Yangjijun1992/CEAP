# μ子 S1 后 Afterpulse 随机耦合本底分析报告

> 场景：RELICS 类探测器，μ 子在液氙中沉积 50 MeV，100 ns 内产生 391,000 PE
> 的快闪 S1 信号；56 只 PMT 均匀接收；全部 PMT APP = 5%。
> 分析窗口：S1 后 [100 ns, 5 μs]；信号符合窗口：1 μs。
> 参数取自 config/settings.yaml（2026-08-17 更新，CE=46%，PDF 50 ns bin）。
> 日期：2026-08-17

---

## 1. 场景参数

| 参数 | 数值 | 说明 |
|------|------|------|
| 液氙能损 | 2 MeV/cm | dE/dx |
| 探测器内径迹长度 | 25 cm | |
| 探测器内沉积能量 | 50 MeV | 25 cm × 2 MeV/cm |
| S1 快闪时长 | 100 ns | |
| 光产额（LXe 本征） | 50 ph/keV | 每 keV 光子数 |
| 量子效率 × 收集效率 | 34% × 46% | 有效探测效率 15.6% |
| 有效光产额 | 7.82 PE/keV | 50 × 0.34 × 0.46 |
| S1 总 PE | 391,000 PE | 50,000 keV × 7.82 PE/keV |
| PMT 数量 | 56 只 | 均匀接收 S1 |
| 每 PMT 主脉冲 PE | ~6,982 PE | 391k / 56 |
| PMT APP | 5% | 全部 PMT，用 LV2358 PDF (app≈5.3%) |
| PDF 时间 bin | 50 ns | 0–60 μs，1200 bin |
| afterpulse 分析窗口 | 100 ns – 5 μs | S1 后 |
| 信号符合窗口 | 1 μs | |
| CEvNS 信号阈值 | 120 – 240 PE | 需求 S-02 |

**模拟工具**：`src/ceap/tools/muon_s1_afterpulse_analysis.py`（配置驱动）
**输出**：`output/muon_s1/muon_s1_afterpulse_analysis.png`、`muon_s1_afterpulse_summary.json`
**统计量**：200 个 μ 子

---

## 2. Afterpulse 电荷与时间分布特征

| 特征 | 数值 | 说明 |
|------|------|------|
| 单 μ 子 afterpulse 事件数 | **8,328 ± 91** | 解析期望 8,543 |
| 单 μ 子 afterpulse 总电荷 | **20,338 ± 344 PE** | 解析期望 ~20,340，偏差 <0.1% |
| 事件电荷均值 | 2.44 PE | 高度偏态分布 |
| 事件电荷中位数 | 1.0 PE | 绝大多数为 1 PE |
| 事件电荷 p99 | 15.0 PE | 长尾至 30 PE (max_npe 上限) |
| 电荷-时间峰值 | Δt ≈ 250 ns | 快衰减 + 离子峰叠加 |
| 前 1 μs 电荷占比 | ~57% | 窗口积分 |
| 前 2 μs 电荷占比 | ~88% | 窗口积分 |

### 电荷分布（事件级）

```
概率密度 (log-log):
  1 PE: 峰值 (中位数)
  均值 2.44 PE
  p99 = 15 PE
  上限 30 PE (max_npe, AP-04)
```

### 时间分布（电荷加权）

```
电荷/bin
    │  ██
    │  ██
    │  ██
    │  ██  ██
    │  ██  ██
    │  ██  ██  ██
    │  ██  ██  ██  ██  ██
    └──────────────────────► Δt (ns)
       0.1  1   2   3   4  5 μs
   快衰减主导，叠加离子峰 (H⁺/He⁺/CH₄⁺ 等)
```

---

## 3. 1 μs 窗口内随机符合事例率

扫描方式：1 μs 滑动窗口，50 ns 步进，覆盖 [100 ns, 5 μs]，每 μ 子 78 个窗口。

| 指标 | 数值 |
|------|------|
| 1 μs 窗内平均 PE | **2,936 PE** |
| 1 μs 窗内 p99 PE | 15,086 PE |
| 窗内总 PE ≥ 120 PE 占比 | **82.73%** |
| 窗内总 PE ∈ [120, 240] 占比（CEvNS 阈值窗命中） | **12.94%** |
| **每 μ 子假信号数**（1 μs 窗∈[120,240] PE） | **10.10 ± 2.48 个** |
| **假信号率**（× μ 子率 R_mu） | **R_mu × 10.10 Hz** |

### 本底率换算示例

| μ 子率 R_mu | 假信号率 |
|-------------|----------|
| 0.01 Hz（100 s 一个 μ 子） | 0.101 Hz |
| 0.1 Hz（10 s 一个 μ 子） | 1.01 Hz |
| 1 Hz（1 s 一个 μ 子） | 10.1 Hz |
| 1.49 Hz（本次配置 47.5×0.0314） | **15.0 Hz**（100k 模拟：14.80 Hz） |

---

## 4. 关键物理结论

1. **afterpulse 电荷巨大**：391k PE 主脉冲 × 5% APP 在 5 μs 内产生约 **2 万 PE**
   afterpulse 电荷，不可忽略。
2. **时间高度集中**：约 88% 的 afterpulse 电荷落在前 2 μs，S1 后紧邻的 1 μs 窗
   平均约 2,900 PE —— **远超任何信号阈值**（CEvNS 120–240 PE）。
3. **假信号率非线性放大**：每 μ 子产生约 **10.1 个**落入 CEvNS 阈值窗的假信号；
   按本次配置 μ 子率 1.49 Hz 折合，假信号率约 **15.0 Hz**（100k 模拟为 14.80 Hz，
   统计涨落范围内一致）。
4. **死区建议**：死区排除必须覆盖 S1 后至少 **2–5 μs**（框架
   `dead_zone.s1_after_ns` 可配），否则 afterpulse 本底会完全淹没 CEvNS 窗口。

---

## 5. 模拟实现与复用

分析脚本：`src/ceap/tools/muon_s1_afterpulse_analysis.py`（配置驱动）

```
场景参数 (来自 config/settings.yaml):
  detector.n_pmts = 56
  detector.single_pe = {quantum_efficiency: 34%, collection_efficiency: 46%}
  detector.target.light_yield_pe_keV = 50 ph/keV
  muon.track_length_mean_cm = 25 cm, dedx = 2 MeV/cm
  signal_window = {length_us: 1.0, threshold_pe_min: 120, threshold_pe_max: 240}
  有效光产额 = 50 × 0.34 × 0.46 = 7.82 PE/keV → S1 = 391,000 PE
  分析窗口: [100 ns, 5 μs]; 统计 μ 子数: 200
```

PDF 数据：`data/app_pdf_LV2358.npz`（APP≈5.3%，代表 5% 情形，时间 bin 50 ns）

输出文件：
- `output/muon_s1/muon_s1_afterpulse_analysis.png`（4 面板：电荷分布 / 时间分布 / 1μs窗PE分布 / 窗内电荷vs位置）
- `output/muon_s1/muon_s1_afterpulse_summary.json`（结构化数值）

---

## 6. 后续可扩展项

| 事项 | 说明 |
|------|------|
| APP 对比扫描 | 5% / 10%（LV2358/LV2305）对比 |
| μ 子率折合绝对率 | 按探测器面积与深度折合 R_mu，给出绝对本底率 (Hz) |
| 死区扫描 | `dead_zone.s1_after_ns` 对假信号率的抑制曲线 |
| 固化为框架配置 | `config/muon_s1_scenario.yaml` 端到端复用 |
| 与纯噪声本底合并 | V 组验证：与暗噪声/物理本底叠加的总本底率 |

---

*本文档对应需求：目标 2（大信号 μ 子 S1 过后 afterpulse 随机耦合本底水平），
模拟结果可由 V 组用于验证与系统误差预算。*
