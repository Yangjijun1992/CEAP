"""μ子 S1 后 afterpulse 随机耦合本底专项分析。

场景 (RELICS 类探测器):
  - μ 子在液氙中沉积 120 MeV → 100 ns 内产生 1,000,000 PE 的快闪 S1
  - 56 只 PMT 均匀接收 → 每只 ~17,857 PE
  - 全部 PMT APP = 5% (用 LV2358 PDF, app≈5.3%)
  - afterpulse 分析窗口: S1 后 [100 ns, 5 μs]
  - 信号符合窗口: 1 μs

输出:
  1) afterpulse 电荷大小分布 + 时间分布特征 (图 + 统计)
  2) 1 μs 窗内随机符合事例率 (含 CEvNS 120-240 PE 阈值)
  3) 事件率分解: 噪声事例率 (Hz) 与 每 μ 子本底事例数

用法: python -m ceap.tools.muon_s1_afterpulse_analysis
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

_CJK = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
if _CJK.exists():
    font_manager.fontManager.addfont(str(_CJK))
    plt.rcParams["font.sans-serif"] = ["Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

from ..models.afterpulse_pdf import FileAfterpulsePDF

OUT_DIR = Path("output/muon_s1")
PDF_FILE = "data/app_pdf_LV2358.npz"

N_PMT = 56
MAIN_PE_TOTAL = 1_000_000
T_WIN_MIN_NS = 100.0        # afterpulse 分析窗口下界
T_WIN_MAX_NS = 5_000.0      # afterpulse 分析窗口上界
WINDOW_US = 1.0             # 信号符合窗口
THR_MIN_PE = 120.0          # CEvNS 阈值下界
THR_MAX_PE = 240.0          # CEvNS 阈值上界
N_MUONS = 200               # 统计用 μ 子数


def run_muon_afterpulse(seed: int):
    """模拟一个 μ 子 S1 的全部 afterpulse 击中, 返回 (dt_ns, npe) 数组。"""
    rng = np.random.default_rng(seed)
    pdf = FileAfterpulsePDF(PDF_FILE, max_npe=30, seed=seed)
    per_pmt = MAIN_PE_TOTAL // N_PMT
    hits = []
    for pmt in range(N_PMT):
        seq = pdf.sample(per_pmt)
        if len(seq):
            hits.append(seq)
    if not hits:
        return np.empty((0, 2))
    return np.vstack(hits)


def analyze_muon(dt_ns, npe):
    """统计单个 μ 子 afterpulse 在 [100ns,5us] 内的特征。"""
    mask = (dt_ns >= T_WIN_MIN_NS) & (dt_ns < T_WIN_MAX_NS)
    dt, npe = dt_ns[mask], npe[mask]

    total_charge = float(npe.sum())
    n_events = len(dt)

    # 时间特征
    dt_edges = np.linspace(T_WIN_MIN_NS, T_WIN_MAX_NS, 50)
    charge_hist, _ = np.histogram(dt, weights=npe, bins=dt_edges)
    dt_c = 0.5 * (dt_edges[:-1] + dt_edges[1:])
    peak_dt_ns = dt_c[int(np.argmax(charge_hist))] if charge_hist.sum() else float("nan")

    # 电荷特征
    charge_mean = float(npe.mean()) if n_events else 0.0
    charge_median = float(np.median(npe)) if n_events else 0.0
    charge_p99 = float(np.percentile(npe, 99)) if n_events else 0.0

    return {
        "n_events": n_events,
        "total_charge": total_charge,
        "charge_mean": charge_mean,
        "charge_median": charge_median,
        "charge_p99": charge_p99,
        "peak_dt_ns": peak_dt_ns,
        "dt": dt,
        "npe": npe,
    }


def coincidence_scan(dt_ns, npe, seed: int):
    """1 μs 滑动窗口扫描 [100ns,5us]: 返回每窗口总 PE 与超阈值标记。"""
    rng = np.random.default_rng(seed)
    # 在 [100ns, 5us-1us] 内放置 1μs 窗口
    t_start = T_WIN_MIN_NS
    t_end = T_WIN_MAX_NS - WINDOW_US * 1e3
    step = 50.0  # ns 步进
    n_win = int((t_end - t_start) / step)
    win_pe = np.zeros(n_win)
    for i in range(n_win):
        lo = t_start + i * step
        hi = lo + WINDOW_US * 1e3
        mask = (dt_ns >= lo) & (dt_ns < hi)
        win_pe[i] = npe[mask].sum()
    thr_hit = (win_pe >= THR_MIN_PE) & (win_pe <= THR_MAX_PE)
    return win_pe, thr_hit


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_dt, all_npe, all_win, all_hit = [], [], [], []
    stats = []
    for i in range(N_MUONS):
        arr = run_muon_afterpulse(seed=1000 + i)
        dt, npe = arr[:, 0], arr[:, 1]
        s = analyze_muon(dt, npe)
        stats.append(s)
        all_dt.append(s["dt"]); all_npe.append(s["npe"])
        win, hit = coincidence_scan(s["dt"], s["npe"], seed=2000 + i)
        all_win.append(win); all_hit.append(hit)

    dt_all = np.concatenate(all_dt)
    npe_all = np.concatenate(all_npe)
    win_all = np.concatenate(all_win)
    hit_all = np.concatenate(all_hit)

    # ---- 统计汇总 ----
    print("=" * 70)
    print("μ子 S1 (1M PE/56 PMT, APP=5%) afterpulse 在 [100ns, 5us] 内特征")
    print("=" * 70)
    print(f"  单 μ 子 afterpulse 事件数    : {np.mean([s['n_events'] for s in stats]):.0f} ± {np.std([s['n_events'] for s in stats]):.0f}")
    print(f"  单 μ 子 afterpulse 总电荷(PE): {np.mean([s['total_charge'] for s in stats]):.0f} ± {np.std([s['total_charge'] for s in stats]):.0f}")
    print(f"  电荷分布: 均值={np.mean(npe_all):.2f} PE  中位数={np.median(npe_all):.2f} PE  p99={np.percentile(npe_all,99):.1f} PE")
    print(f"  峰值时间(dt)               : {np.median([s['peak_dt_ns'] for s in stats]):.0f} ns")

    print(f"\n  --- 1 μs 窗内随机符合 (步进 50 ns, 每 μ 子扫描窗数 {win_all.size//N_MUONS:.0f}) ---")
    print(f"  全部窗平均 PE     : {win_all.mean():.1f} PE")
    print(f"  全部窗 99% 分位   : {np.percentile(win_all,99):.1f} PE")
    frac_120 = (win_all >= THR_MIN_PE).mean()
    frac_120_240 = hit_all.mean()
    print(f"  窗内总PE≥120 PE   : {frac_120*100:.2f}%")
    print(f"  窗内总PE∈[120,240]: {frac_120_240*100:.3f}%  (CEvNS 阈值窗口命中率)")

    # 每 μ 子命中数 (一个 μ 子的 afterpulse 产生的假信号数)
    per_muon = np.array([hit.sum() for hit in np.split(hit_all, N_MUONS)])
    print(f"  每 μ 子假信号数   : {per_muon.mean():.2f} ± {per_muon.std():.2f} (1μs窗∈[120,240]PE)")

    # 换算事例率: 假设 μ 子率 R_mu, 假信号率 = R_mu × 每μ子假信号数
    # (用户未给 μ 子率, 给出归一化结果 + 与窗内 PE 的对应)
    print(f"\n  --- 本底率换算 (× μ子率 R_mu) ---")
    print(f"  假信号率 = R_mu × {per_muon.mean():.3f} Hz/(μ子/s)")
    print(f"  若 R_mu = 0.01 Hz (100 s 一个 μ子): {0.01*per_muon.mean():.4f} Hz")

    # ---- 绘图 ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle("μ子 S1 (1M PE/56 PMT, APP=5%) afterpulse 特征", fontsize=14)

    # 1. 电荷分布
    ax = axes[0, 0]
    bins = np.arange(0.5, 20.5, 1.0)
    ax.hist(npe_all, bins=bins, density=True, alpha=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("afterpulse 电荷 (PE)")
    ax.set_ylabel("概率密度")
    ax.set_title(f"电荷分布 (N={npe_all.size:.0f}, 均值 {np.mean(npe_all):.2f} PE)")

    # 2. 时间分布 (电荷加权)
    ax = axes[0, 1]
    dt_edges = np.linspace(T_WIN_MIN_NS, T_WIN_MAX_NS, 100)
    ax.hist(dt_all, bins=dt_edges, weights=npe_all, alpha=0.8)
    ax.set_xlabel("Δt (ns after S1)")
    ax.set_ylabel("电荷 (PE/bin)")
    ax.set_title("afterpulse 电荷-时间分布")

    # 3. 1μs 窗内 PE 分布
    ax = axes[1, 0]
    bins = np.logspace(np.log10(max(win_all.min(),1)), np.log10(win_all.max()), 40)
    ax.hist(win_all, bins=bins, alpha=0.8)
    ax.axvline(THR_MIN_PE, color="r", ls="--", label=f"阈值 {THR_MIN_PE:.0f} PE")
    ax.axvline(THR_MAX_PE, color="r", ls=":", label=f"阈值 {THR_MAX_PE:.0f} PE")
    ax.set_xscale("log")
    ax.set_xlabel("1 μs 窗内总 PE")
    ax.set_ylabel("窗数")
    ax.set_title("1 μs 窗口内随机符合 PE 分布")
    ax.legend()

    # 4. 窗内 PE vs 窗口位置 (单 μ 子示例)
    ax = axes[1, 1]
    dt_edges = np.linspace(T_WIN_MIN_NS, T_WIN_MAX_NS, 100)
    charge_hist, _ = np.histogram(dt_all, weights=npe_all, bins=dt_edges)
    dt_c = 0.5 * (dt_edges[:-1] + dt_edges[1:])
    # 1μs 窗内电荷 (滑动平均近似)
    win_charge = np.convolve(charge_hist, np.ones(20), mode="same") / 20
    ax.plot(dt_c, win_charge, lw=1.5)
    ax.axhline(THR_MIN_PE, color="r", ls="--", label=f"阈值 {THR_MIN_PE:.0f} PE")
    ax.axhline(THR_MAX_PE, color="r", ls=":", label=f"阈值 {THR_MAX_PE:.0f} PE")
    ax.set_xlabel("Δt (ns after S1)")
    ax.set_ylabel("1 μs 窗内电荷 (PE)")
    ax.set_yscale("log")
    ax.set_title("1 μs 窗内电荷 vs 窗口位置")
    ax.legend()

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_DIR / "muon_s1_afterpulse_analysis.png", dpi=130)
    plt.close(fig)

    # ---- JSON ----
    result = {
        "scenario": {
            "n_pmt": N_PMT, "main_pe_total": MAIN_PE_TOTAL,
            "main_pe_per_pmt": MAIN_PE_TOTAL // N_PMT,
            "app": 0.053, "pdf": PDF_FILE,
            "analysis_window_ns": [T_WIN_MIN_NS, T_WIN_MAX_NS],
            "coincidence_window_us": WINDOW_US,
        },
        "per_muon": {
            "n_events_mean": float(np.mean([s["n_events"] for s in stats])),
            "n_events_std": float(np.std([s["n_events"] for s in stats])),
            "total_charge_mean_pe": float(np.mean([s["total_charge"] for s in stats])),
            "total_charge_std_pe": float(np.std([s["total_charge"] for s in stats])),
            "charge_mean_pe": float(np.mean(npe_all)),
            "charge_median_pe": float(np.median(npe_all)),
            "charge_p99_pe": float(np.percentile(npe_all, 99)),
            "peak_dt_ns": float(np.median([s["peak_dt_ns"] for s in stats])),
        },
        "coincidence_1us": {
            "window_mean_pe": float(win_all.mean()),
            "window_p99_pe": float(np.percentile(win_all, 99)),
            "frac_ge_120pe": float(frac_120),
            "frac_120_240pe": float(frac_120_240),
            "fake_per_muon": float(per_muon.mean()),
            "fake_per_muon_std": float(per_muon.std()),
        },
    }
    with open(OUT_DIR / "muon_s1_afterpulse_summary.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n结果已保存: {OUT_DIR}/")
    return result


if __name__ == "__main__":
    main()
