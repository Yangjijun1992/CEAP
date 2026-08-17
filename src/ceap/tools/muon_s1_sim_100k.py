"""10 万次 μ子 S1 afterpulse 随机耦合本底大规模模拟 (向量化批处理)。

场景 (RELICS, 参数取自 config/settings.yaml):
  - μ子 S1: 由 muon 模型计算 (能损 × 光产额 × QE×CE) / 56 PMT
  - μ子在探测器中的事例率 = 通量 × 有效面积 (可由 CLI 参数覆盖)
  - afterpulse 分析窗口: S1 后 [100 ns, 5 μs]
  - 信号符合窗口: 1 μs (滑动, 步进 50 ns)
  - 假信号定义: 1 μs 窗内总 PE ∈ [threshold_pe_min, threshold_pe_max]

实现:
  - 每 batch 500 μ子, 二项抽样事件数, 向量化从 2D PDF 抽样 (dt, npe)
  - 按 50 ns bin 累加电荷谱, 滑动窗口求和, 判定假信号
  - 统计: 假信号 PE / 时间宽度 / 距 μ子 Δt

输出 (output/muon_s1_100k/):
  - fake_signal_energy_spectrum.png  完整假信号能谱 (log-log, 全部窗 + 120/240 边界)
  - fake_signal_energy_compare.png   线性轴对比: 所有窗 vs 假信号 dN/dPE + 幂律拟合
  - fake_signal_pe.png       信号窗内假信号 PE 直方图 (阈值红线)
  - fake_signal_width.png    假信号时间宽度直方图
  - fake_signal_delta_t.png  假信号距 μ子 Δt 直方图 (μs)
  - summary.json

用法: python -m ceap.tools.muon_s1_sim_100k [--n-muons 100000] [--rate HZ] [--config FILE]
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

from ..config.loader import load_config
from ..models.afterpulse_pdf import FileAfterpulsePDF
from ..models.muon import MuonModel

DEFAULT_PDF_FILE = "data/app_pdf_LV2358.npz"   # APP≈5% 参考 PDF
T_MIN_NS = 100.0
T_MAX_NS = 5_000.0
STEP_NS = 50.0
BIN_NS = 50.0
BATCH = 500


def load_scenario(cfg_path: str | None) -> dict:
    """从配置派生场景参数。"""
    cfg = load_config(cfg_path)
    det = cfg.get("detector") or {}
    sw = cfg.get("signal_window") or {}
    ap = cfg.get("afterpulse") or {}
    muon = MuonModel(cfg)

    n_pmt = int(det.get("n_pmts", 0)) or 56
    main_pe_total = int(round(muon.mean_energy_mev() * 1e3 * muon.ly))
    pdf_file = ap.get("pdf_file") or DEFAULT_PDF_FILE

    return {
        "n_pmt": n_pmt,
        "main_pe_total": main_pe_total,
        "main_pe_per_pmt": main_pe_total // n_pmt,
        "muon_rate_hz": muon.hit_rate_hz,
        "pdf_file": pdf_file,
        "window_ns": float(sw.get("length_us", 1.0)) * 1e3,
        "thr_min": float(sw.get("threshold_pe_min", 120.0)),
        "thr_max": float(sw.get("threshold_pe_max", 240.0)),
    }


def _sample_batch(pdf: FileAfterpulsePDF, rng: np.random.Generator, n_mu: int, n_pmt: int, per_pmt: int):
    """一批 n_mu 个 μ子的 (muon_id, dt_ns, npe)。"""
    n_ev = rng.binomial(per_pmt, pdf.p_ap_total, size=(n_mu, n_pmt))
    tot = int(n_ev.sum())
    muon_id = np.repeat(np.arange(n_mu), n_ev.sum(axis=1))
    pmt_id = np.repeat(np.tile(np.arange(n_pmt), n_mu), n_ev.ravel())
    idx = rng.choice(pdf._flat_idx, size=tot, p=pdf._flat_pdf[pdf._flat_idx])
    i_dt, i_npe = np.unravel_index(idx, pdf.hist.shape)
    dt = pdf.dt_centers[i_dt]
    npe = np.minimum(pdf.npe_centers[i_npe], float(pdf.max_npe))
    return muon_id, pmt_id, dt, npe


def simulate(sc: dict, n_muons: int, seed: int = 42):
    pdf = FileAfterpulsePDF(sc["pdf_file"], max_npe=30, seed=seed)
    rng = np.random.default_rng(seed)
    n_pmt, per_pmt = sc["n_pmt"], sc["main_pe_per_pmt"]
    window_ns, thr_min, thr_max = sc["window_ns"], sc["thr_min"], sc["thr_max"]

    n_win = int((T_MAX_NS - T_MIN_NS - window_ns) / STEP_NS) + 1
    win_start = T_MIN_NS + np.arange(n_win) * STEP_NS
    span = int(window_ns / BIN_NS)
    n_bins = int((T_MAX_NS - T_MIN_NS) / BIN_NS)

    pe_list, width_list, dt_list = [], [], []
    win_all_list = []
    n_fake_total = 0
    n_afterpulse_total = 0

    t0 = time.time()
    for start in range(0, n_muons, BATCH):
        n_mu = min(BATCH, n_muons - start)
        muon_id, pmt_id, dt, npe = _sample_batch(pdf, rng, n_mu, n_pmt, per_pmt)
        n_afterpulse_total += len(dt)

        # 50ns bin 累加电荷谱 (n_mu, n_bins)
        bin_idx = ((dt - T_MIN_NS) / BIN_NS).astype(int)
        ok = (bin_idx >= 0) & (bin_idx < n_bins)
        muon_id, bin_idx, npe = muon_id[ok], bin_idx[ok], npe[ok]
        flat = muon_id * n_bins + bin_idx
        charge = np.bincount(flat, weights=npe, minlength=n_mu * n_bins).reshape(n_mu, n_bins)

        # 滑动窗口: 窗 i 覆盖 bins [i, i+span-1], window_ns 窗和 = span 个 bin 之和
        win_sum = np.zeros((n_mu, n_win))
        for i in range(n_win):
            win_sum[:, i] = charge[:, i:i + span].sum(axis=1)
        win_all_list.append(win_sum.astype(np.float32))

        mask = (win_sum >= thr_min) & (win_sum <= thr_max)
        idx_mu, idx_w = np.nonzero(mask)
        if len(idx_mu):
            pe = win_sum[idx_mu, idx_w]
            span_mask = charge[idx_mu][:, :span] > 0
            first = np.argmax(span_mask, axis=1)
            last = span - 1 - np.argmax(span_mask[:, ::-1], axis=1)
            width_ns = (last - first + 1) * BIN_NS
            dt_win = win_start[idx_w] + window_ns / 2
            pe_list.append(pe); width_list.append(width_ns); dt_list.append(dt_win)
            n_fake_total += len(pe)

    pe = np.concatenate(pe_list) if pe_list else np.empty(0)
    width = np.concatenate(width_list) if width_list else np.empty(0)
    dt_win = np.concatenate(dt_list) if dt_list else np.empty(0)
    win_all = np.concatenate(win_all_list) if win_all_list else np.empty(0)
    win_all = np.asarray(win_all).ravel()

    # 事例率: 模拟挂机时间 = n_muons / rate_hz
    live_time_s = n_muons / sc["muon_rate_hz"]
    rate_hz = n_fake_total / live_time_s

    return {
        "n_muons": n_muons,
        "n_afterpulse_total": n_afterpulse_total,
        "n_fake": n_fake_total,
        "fake_per_muon": n_fake_total / n_muons,
        "rate_hz": rate_hz,
        "live_time_s": live_time_s,
        "pe": pe, "width_ns": width, "dt_ns": dt_win,
        "win_all": win_all,
        "elapsed_s": time.time() - t0,
    }


def plot_histograms(res: dict, sc: dict, out_dir: Path):
    OUT_DIR = Path(out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    thr_min, thr_max = sc["thr_min"], sc["thr_max"]
    window_ns = sc["window_ns"]

    # 1. 假信号 PE 直方图 (阈值红线, 边框标注)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["pe"], bins=60, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(thr_min, color="red", ls="--", lw=1.8)
    ax.axvline(thr_max, color="red", ls="--", lw=1.8)
    ax.axvspan(thr_min, thr_max, color="red", alpha=0.12)
    ax.annotate(f"{thr_min:.0f} PE（下界）", xy=(thr_min, ax.get_ylim()[1]),
                xytext=(thr_min + 0.01 * (thr_max - thr_min), ax.get_ylim()[1] * 0.97),
                color="red", fontsize=10, ha="left", va="top")
    ax.annotate(f"{thr_max:.0f} PE（上界）", xy=(thr_max, ax.get_ylim()[1]),
                xytext=(thr_max - 0.01 * (thr_max - thr_min), ax.get_ylim()[1] * 0.97),
                color="red", fontsize=10, ha="right", va="top")
    ax.set_xlabel(f"假信号总 PE / 能量 ({window_ns/1e3:.0f} μs 窗)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title(f"假信号能量谱 (N={len(res['pe']):,}, 均值 {res['pe'].mean():.1f} PE)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_pe.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 2. 假信号时间宽度直方图
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["width_ns"] / 1e3, bins=40, color="seagreen", edgecolor="white", alpha=0.85)
    ax.set_xlabel("假信号时间宽度 (μs)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title(f"假信号时间宽度分布 (均值 {res['width_ns'].mean()/1e3:.3f} μs)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_width.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 3. 假信号距 μ子 Δt 直方图 (μs)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.grid(True, alpha=0.4)
    ax.hist(res["dt_ns"] / 1e3, bins=80, color="darkorange", edgecolor="white", alpha=0.85)
    ax.set_xlabel("距 μ子 S1 的 Δt (μs)", fontsize=11)
    ax.set_ylabel("计数", fontsize=11)
    ax.set_title("假信号 Δt 分布 (窗中心)", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_delta_t.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_full_energy_spectrum(res: dict, sc: dict, out_dir: Path):
    """完整假信号能谱: 全部 1 μs 窗总 PE 分布 (log-log), 标注 120/240 边界。"""
    OUT_DIR = Path(out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    thr_min, thr_max = sc["thr_min"], sc["thr_max"]
    window_ns = sc["window_ns"]
    win_all = res["win_all"]

    n_total = win_all.size
    n_in_band = int(((win_all >= thr_min) & (win_all <= thr_max)).sum())
    n_above = int((win_all > thr_max).sum())

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.grid(True, alpha=0.4, which="both")
    # 全部窗电荷分布 (log-log)
    lo = max(float(win_all.min()), 1.0)
    hi = float(win_all.max())
    bins = np.logspace(np.log10(lo), np.log10(hi), 80)
    ax.hist(win_all, bins=bins, color="steelblue", edgecolor="white", alpha=0.9,
            label=f"全部窗 (N={n_total:,})")
    # 信号区间
    ax.axvline(thr_min, color="red", ls="--", lw=1.8)
    ax.axvline(thr_max, color="red", ls="--", lw=1.8)
    ax.axvspan(thr_min, thr_max, color="red", alpha=0.15)
    ax.annotate(f"{thr_min:.0f} PE（下界）", xy=(thr_min, ax.get_ylim()[1]),
                xytext=(thr_min * 1.06, ax.get_ylim()[1] * 0.96),
                color="red", fontsize=10, ha="left", va="top")
    ax.annotate(f"{thr_max:.0f} PE（上界）", xy=(thr_max, ax.get_ylim()[1]),
                xytext=(thr_max * 0.8, ax.get_ylim()[1] * 0.88),
                color="red", fontsize=10, ha="center", va="top")
    # 假信号 (落在信号窗内的窗) 叠加
    ax.hist(res["pe"], bins=np.linspace(thr_min, thr_max, 30), color="red",
            alpha=0.55, label=f"假信号 (N={n_in_band:,}, {n_in_band/n_total*100:.2f}%)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(f"1 μs 窗内总 PE / 能量 ({window_ns/1e3:.0f} μs 窗)", fontsize=12)
    ax.set_ylabel("窗口计数", fontsize=12)
    ax.set_title(f"完整假信号能谱 (log-log)\n"
                 f"信号窗 [{thr_min:.0f}, {thr_max:.0f}] PE 内 {n_in_band:,} 窗 ({n_in_band/n_total*100:.2f}%)；"
                 f"> {thr_max:.0f} PE 共 {n_above:,} 窗 ({n_above/n_total*100:.2f}%)",
                 fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_energy_spectrum.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  完整能谱: 总窗 {n_total:,}, 信号窗内 {n_in_band:,} ({n_in_band/n_total*100:.2f}%), "
          f">{thr_max:.0f}PE {n_above:,} ({n_above/n_total*100:.2f}%)")


def plot_linear_comparison(res: dict, sc: dict, out_dir: Path):
    """统一线性轴对比: 所有窗 dN/dPE vs 假信号 dN/dPE + 幂律拟合。"""
    OUT_DIR = Path(out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    thr_min, thr_max = sc["thr_min"], sc["thr_max"]
    win = res["win_all"]
    fake = res["pe"]

    # 统一线性 bin (X), 取 dN/dPE (每 PE 计数)
    bin_w = 20.0
    xmax = 1000.0
    edges = np.arange(0.0, xmax + bin_w, bin_w)
    h_all, _ = np.histogram(win, bins=edges)
    dndpe_all = h_all / bin_w
    h_f, _ = np.histogram(fake, bins=edges)
    dndpe_f = h_f / bin_w
    xc = (edges[:-1] + edges[1:]) / 2.0

    # 幂律拟合: dN/dPE = A * PE^(-alpha), 用 log 密度
    win_all = np.asarray(win).ravel()
    f_edges = np.logspace(np.log10(max(thr_min, 1.0)), np.log10(win_all.max()), 40)
    fh, _ = np.histogram(win_all, bins=f_edges)
    fd = fh / np.diff(f_edges)
    fxc = (f_edges[:-1] + f_edges[1:]) / 2.0
    m = fd > 0
    fit = np.polyfit(np.log(fxc[m]), np.log(fd[m]), 1)
    alpha = -fit[0]
    A = np.exp(fit[1])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.grid(True, alpha=0.4)
    # 所有窗 dN/dPE (曲线)
    ax.plot(xc, dndpe_all, "-", color="steelblue", lw=2.2,
            label=f"所有窗 dN/dPE (N={win.size:,})")
    # 假信号 dN/dPE (同一线性轴, 应落在 [120,240])
    ax.plot(xc, dndpe_f, "-o", color="crimson", ms=3, lw=1.6,
            label=f"假信号 dN/dPE (N={len(fake):,}, ∈[{thr_min:.0f},{thr_max:.0f}])")
    # 幂律拟合线 (在拟合坐标上画)
    fxs = np.logspace(np.log10(max(thr_min, 1.0)), np.log10(win_all.max()), 100)
    ax.plot(fxs, A * np.power(fxs, -alpha), "--", color="green", lw=1.8,
            label=f"幂律拟合 dN/dPE ≈ {A:.2g}·PE^({-alpha:.2f})")
    # 信号窗边框
    ax.axvspan(thr_min, thr_max, color="red", alpha=0.12)
    ax.axvline(thr_min, color="red", ls="--", lw=1.2)
    ax.axvline(thr_max, color="red", ls="--", lw=1.2)
    ax.annotate(f"{thr_min:.0f} PE", xy=(thr_min, ax.get_ylim()[1]),
                xytext=(thr_min, ax.get_ylim()[1] * 0.98), color="red",
                fontsize=10, ha="left", va="top")
    ax.annotate(f"{thr_max:.0f} PE", xy=(thr_max, ax.get_ylim()[1]),
                xytext=(thr_max, ax.get_ylim()[1] * 0.98), color="red",
                fontsize=10, ha="right", va="top")
    ax.set_xlim(0, xmax)
    ax.set_xlabel("1 μs 窗内总 PE（线性轴）", fontsize=12)
    ax.set_ylabel("dN/dPE（每 PE 窗口计数，线性轴）", fontsize=12)
    ax.set_title(f"所有窗 vs 假信号能谱对比（统一线性轴）\n"
                 f"幂指数 α={alpha:.3f}：dN/dPE ∝ PE^(-{alpha:.2f})", fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fake_signal_energy_compare.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"线性对比图已保存; 幂指数 dN/dPE ∝ PE^(-{alpha:.3f}) (A={A:.3g})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-muons", type=int, default=100_000)
    ap.add_argument("--rate", type=float, default=None,
                    help="μ子事例率 (默认取 config 通量×面积)")
    ap.add_argument("--config", default=None, help="配置文件 (默认 config/settings.yaml)")
    args = ap.parse_args()

    sc = load_scenario(args.config)
    if args.rate:
        sc["muon_rate_hz"] = args.rate

    res = simulate(sc, args.n_muons)
    out = Path("output/muon_s1_100k")
    plot_histograms(res, sc, out)
    plot_full_energy_spectrum(res, sc, out)
    plot_linear_comparison(res, sc, out)

    summary = {
        "scenario": sc,
        "n_muons": res["n_muons"],
        "muon_rate_hz": sc["muon_rate_hz"],
        "live_time_s": res["n_muons"] / sc["muon_rate_hz"],
        "afterpulse_events_total": res["n_afterpulse_total"],
        "afterpulse_per_muon": res["n_afterpulse_total"] / res["n_muons"],
        "fake_signals_total": res["n_fake"],
        "fake_per_muon": res["fake_per_muon"],
        "fake_rate_hz": res["n_fake"] / (res["n_muons"] / sc["muon_rate_hz"]),
        "fake_pe_mean": float(res["pe"].mean()) if len(res["pe"]) else 0.0,
        "fake_pe_median": float(np.median(res["pe"])) if len(res["pe"]) else 0.0,
        "fake_width_mean_us": float(res["width_ns"].mean() / 1e3) if len(res["width_ns"]) else 0.0,
        "fake_dt_mean_us": float(res["dt_ns"].mean() / 1e3) if len(res["dt_ns"]) else 0.0,
        "window_stats": {
            "n_windows_total": int(res["win_all"].size),
            "n_in_signal_window": int(((res["win_all"] >= sc["thr_min"]) & (res["win_all"] <= sc["thr_max"])).sum()),
            "frac_in_signal_window": float(((res["win_all"] >= sc["thr_min"]) & (res["win_all"] <= sc["thr_max"])).mean()),
            "n_above": int((res["win_all"] > sc["thr_max"]).sum()),
            "frac_above": float((res["win_all"] > sc["thr_max"]).mean()),
            "n_below": int((res["win_all"] < sc["thr_min"]).sum()),
            "frac_below": float((res["win_all"] < sc["thr_min"]).mean()),
        },
        "elapsed_s": res["elapsed_s"],
    }
    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"plots saved to {out}/")
    return summary


if __name__ == "__main__":
    main()
