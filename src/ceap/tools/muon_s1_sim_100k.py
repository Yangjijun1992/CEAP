"""10 万次 μ子 S1 afterpulse 随机耦合本底大规模模拟 (向量化批处理)。

场景 (RELICS):
  - μ子 S1: 1M PE / 56 PMT (每只 17,857 PE), 全部 PMT APP=5% (LV2358)
  - μ子在探测器中的事例率: 10 Hz (可由 CLI 参数覆盖)
  - afterpulse 分析窗口: S1 后 [100 ns, 5 μs]
  - 信号符合窗口: 1 μs (滑动, 步进 50 ns)
  - 假信号定义: 1 μs 窗内总 PE ∈ [120, 240] (CEvNS 阈值窗)

实现:
  - 每 batch 500 μ子, 二项抽样事件数, 向量化从 2D PDF 抽样 (dt, npe)
  - 按 50 ns bin 累加电荷谱, 滑动窗口求和, 判定假信号
  - 统计: 假信号 PE / 时间宽度 / 距 μ子 Δt

输出 (output/muon_s1_100k/):
  - fake_signal_pe.png       假信号 PE 直方图 (120/240 红线)
  - fake_signal_width.png    假信号时间宽度直方图
  - fake_signal_delta_t.png  假信号距 μ子 Δt 直方图 (μs)
  - summary.json

用法: python -m ceap.tools.muon_s1_sim_100k [--n-muons 100000] [--rate 10]
"""
from __future__ import annotations

import argparse
import json
import time
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

PDF_FILE = "data/app_pdf_LV2358.npz"   # APP≈5%
N_PMT = 56
MAIN_PE_TOTAL = 1_000_000
T_MIN_NS = 100.0
T_MAX_NS = 5_000.0
WINDOW_NS = 1_000.0
STEP_NS = 50.0
THR_MIN = 120.0
THR_MAX = 240.0
BIN_NS = 50.0
BATCH = 500


def _sample_batch(pdf: FileAfterpulsePDF, rng: np.random.Generator, n_mu: int):
    """一批 n_mu 个 μ子的 (muon_id, dt_ns, npe)。"""
    per_pmt = MAIN_PE_TOTAL // N_PMT
    n_ev = rng.binomial(per_pmt, pdf.p_ap_total, size=(n_mu, N_PMT))
    tot = int(n_ev.sum())
    muon_id = np.repeat(np.arange(n_mu), n_ev.sum(axis=1))
    pmt_id = np.repeat(np.tile(np.arange(N_PMT), n_mu), n_ev.ravel())
    idx = rng.choice(pdf._flat_idx, size=tot, p=pdf._flat_pdf[pdf._flat_idx])
    i_dt, i_npe = np.unravel_index(idx, pdf.hist.shape)
    dt = pdf.dt_centers[i_dt]
    npe = np.minimum(pdf.npe_centers[i_npe], float(pdf.max_npe))
    return muon_id, pmt_id, dt, npe


def simulate(n_muons: int, seed: int = 42):
    pdf = FileAfterpulsePDF(PDF_FILE, max_npe=30, seed=seed)
    rng = np.random.default_rng(seed)

    n_win = int((T_MAX_NS - T_MIN_NS - WINDOW_NS) / STEP_NS) + 1
    win_start = T_MIN_NS + np.arange(n_win) * STEP_NS
    span = int(WINDOW_NS / BIN_NS)
    n_bins = int((T_MAX_NS - T_MIN_NS) / BIN_NS)

    pe_list, width_list, dt_list = [], [], []
    n_fake_total = 0
    n_afterpulse_total = 0

    t0 = time.time()
    for start in range(0, n_muons, BATCH):
        n_mu = min(BATCH, n_muons - start)
        muon_id, pmt_id, dt, npe = _sample_batch(pdf, rng, n_mu)
        n_afterpulse_total += len(dt)

        # 50ns bin 累加电荷谱 (n_mu, n_bins)
        bin_idx = ((dt - T_MIN_NS) / BIN_NS).astype(int)
        ok = (bin_idx >= 0) & (bin_idx < n_bins)
        muon_id, bin_idx, npe = muon_id[ok], bin_idx[ok], npe[ok]
        flat = muon_id * n_bins + bin_idx
        charge = np.bincount(flat, weights=npe, minlength=n_mu * n_bins).reshape(n_mu, n_bins)

        # 滑动窗口: 窗 i 覆盖 bins [i, i+span-1], 1μs 窗和 = span 个 bin 之和
        win_sum = np.zeros((n_mu, n_win))
        for i in range(n_win):
            win_sum[:, i] = charge[:, i:i + span].sum(axis=1)

        mask = (win_sum >= THR_MIN) & (win_sum <= THR_MAX)
        idx_mu, idx_w = np.nonzero(mask)
        if len(idx_mu):
            pe = win_sum[idx_mu, idx_w]
            span_mask = charge[idx_mu][:, :span] > 0
            first = np.argmax(span_mask, axis=1)
            last = span - 1 - np.argmax(span_mask[:, ::-1], axis=1)
            width_ns = (last - first + 1) * BIN_NS
            dt_win = win_start[idx_w] + WINDOW_NS / 2
            pe_list.append(pe); width_list.append(width_ns); dt_list.append(dt_win)
            n_fake_total += len(pe)

    pe = np.concatenate(pe_list) if pe_list else np.empty(0)
    width = np.concatenate(width_list) if width_list else np.empty(0)
    dt_win = np.concatenate(dt_list) if dt_list else np.empty(0)

    # 事例率: 模拟挂机时间 = n_muons / rate_hz
    live_time_s = n_muons / 10.0
    rate_hz = n_fake_total / live_time_s

    return {
        "n_muons": n_muons,
        "n_afterpulse_total": n_afterpulse_total,
        "n_fake": n_fake_total,
        "fake_per_muon": n_fake_total / n_muons,
        "rate_hz": rate_hz,
        "live_time_s": live_time_s,
        "pe": pe, "width_ns": width, "dt_ns": dt_win,
        "elapsed_s": time.time() - t0,
    }


def plot_histograms(res: dict, out_dir: Path):
    OUT_DIR = Path(out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. 假信号 PE 直方图 (120/240 红线)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["pe"], bins=np.arange(100, 261, 2), color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(THR_MIN, color="red", ls="--", lw=1.8, label=f"信号区间 {THR_MIN:.0f}–{THR_MAX:.0f} PE")
    ax.axvline(THR_MAX, color="red", ls="--", lw=1.8)
    ax.axvspan(THR_MIN, THR_MAX, color="red", alpha=0.12)
    ax.set_xlabel("假信号总 PE (1 μs 窗)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title(f"假信号 PE 分布 (N={len(res['pe']):,})", fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_pe.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 2. 假信号时间宽度直方图
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["width_ns"] / 1e3, bins=np.arange(0, 1.05, 0.02), color="seagreen", edgecolor="white", alpha=0.85)
    ax.set_xlabel("假信号时间宽度 (μs)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title(f"假信号时间宽度分布 (均值 {res['width_ns'].mean()/1e3:.3f} μs)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_width.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 3. 假信号距 μ子 Δt 直方图 (μs)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["dt_ns"] / 1e3, bins=np.arange(0.1, 4.6, 0.05), color="darkorange", edgecolor="white", alpha=0.85)
    ax.set_xlabel("距 μ子 S1 的 Δt (μs)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title("假信号 Δt 分布 (1 μs 窗中心)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_delta_t.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-muons", type=int, default=100_000)
    ap.add_argument("--rate", type=float, default=10.0)
    args = ap.parse_args()

    res = simulate(args.n_muons)
    out = Path("output/muon_s1_100k")
    plot_histograms(res, out)

    summary = {
        "n_muons": res["n_muons"],
        "muon_rate_hz": args.rate,
        "live_time_s": res["n_muons"] / args.rate,
        "afterpulse_events_total": res["n_afterpulse_total"],
        "afterpulse_per_muon": res["n_afterpulse_total"] / res["n_muons"],
        "fake_signals_total": res["n_fake"],
        "fake_per_muon": res["fake_per_muon"],
        "fake_rate_hz": res["n_fake"] / (res["n_muons"] / args.rate),
        "fake_pe_mean": float(res["pe"].mean()) if len(res["pe"]) else 0.0,
        "fake_pe_median": float(np.median(res["pe"])) if len(res["pe"]) else 0.0,
        "fake_width_mean_us": float(res["width_ns"].mean() / 1e3) if len(res["width_ns"]) else 0.0,
        "fake_dt_mean_us": float(res["dt_ns"].mean() / 1e3) if len(res["dt_ns"]) else 0.0,
        "elapsed_s": res["elapsed_s"],
    }
    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"plots saved to {out}/")
    return summary


if __name__ == "__main__":
    main()
