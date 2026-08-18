"""APP area-vs-delta_time 二维 PDF 可视化 (风格参考 pmt_analysis app.py)。

参考:
  - plot_afterpulse_2d_histogram(): ax1.hist2d(..., bins=[80,80], cmap="jet",
    density=True, norm=LogNorm()) + 下方 Δt 1D 投影 + 离子峰标记
  - plot_afterpulse_delta_time_all_channels(): 3×3 网格 Δt 分布, 时间单位 μs

输出 (output/plots/):
  - app_pdf_<pmt>_2d_before_after.png: 每 PMT 归一化前后对比 (2D + Δt 1D)
  - app_pdf_delta_time_3x3.png: 三个 PMT Δt 分布 3×3 网格 (μs)

用法: python -m ceap.tools.plot_app_pdfs
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import font_manager

_CJK_FONT = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
if _CJK_FONT.exists():
    font_manager.fontManager.addfont(str(_CJK_FONT))
    plt.rcParams["font.sans-serif"] = ["Noto Sans CJK JP", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

from .app_pdf_builder import load_app_npz

RAW_DIR = Path("/mnt/data/PMT/R8520_406/output/app_distributions")
PDF_DIR = Path("data")
OUT_DIR = Path("output/plots")

PMTS = {
    "LV2305": ("app_dist_00381_LV2305.npz", "app_pdf_LV2305.npz", "中 APP ≈ 10%"),
    "LV2358": ("app_dist_00376_LV2358.npz", "app_pdf_LV2358.npz", "低 APP ≈ 5%"),
}


def _scatter(pmt: str):
    """返回原始散点 (delta_time_us, area_pe)。"""
    data = load_app_npz(RAW_DIR / PMTS[pmt][0])
    dt = np.asarray(data["ap_delta_time_ns"], dtype=float) / 1000.0
    area = np.asarray(data["ap_area_pe"], dtype=float)
    mask = np.isfinite(dt) & np.isfinite(area) & (dt >= 0) & (area > 0)
    return dt[mask], area[mask]


def _load_pdf(pmt: str):
    data = np.load(PDF_DIR / PMTS[pmt][1])
    return (
        np.asarray(data["h"], dtype=float),
        np.asarray(data["dt_edges"], dtype=float) / 1000.0,
        np.asarray(data["npe_edges"], dtype=float),
        float(data["p_ap_total"]),
        float(data["p_ap_1us"]),
        float(data["app"]),
    )


def _panel(ax, dt_us, area_pe, density: bool, title: str, dt_max: float, area_max: float):
    """单面板: hist2d + 色标 (风格同 app.py plot_afterpulse_2d_histogram)。"""
    hist = ax.hist2d(
        dt_us, area_pe,
        bins=[80, 80],
        range=[[0, dt_max], [0, area_max]],
        cmap="jet",
        density=density,
        norm=matplotlib.colors.LogNorm(),
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(hist[3], ax=cax)
    cbar.set_label("Density" if density else "Counts", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_ylabel("Afterpulse Area [PE]", fontsize=11)
    ax.tick_params(axis="y", direction="in", labelsize=9, pad=4, length=4, width=1)
    ax.set_title(title, fontsize=11)


def _dt_projection(ax, dt_us, dt_max: float, title: str):
    """下方面板: Δt 1D 直方图 (时间单位 μs)。"""
    ax.hist(
        dt_us,
        bins=80,
        range=(0, dt_max),
        histtype="stepfilled",
        color="skyblue",
        edgecolor="skyblue",
        linewidth=0.8,
        alpha=0.9,
    )
    ax.set_ylabel("Counts", fontsize=11)
    ax.set_xlabel("Time Delay [$\mu$s]", fontsize=11)
    ax.tick_params(axis="x", direction="in", labelsize=9, pad=4, length=4, width=1)
    ax.tick_params(axis="y", direction="in", labelsize=9, pad=4, length=4, width=1)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="y", scilimits=(4, 4))
    ax.yaxis.get_offset_text().set_fontsize(9)
    ax.set_title(title, fontsize=11)


def plot_2d_before_after(pmt: str, label: str):
    """每 PMT: 左列归一化前, 右列归一化后; 每列上下 = 2D + Δt 1D。"""
    dt_us, area_pe = _scatter(pmt)
    h, dt_edges_us, npe_edges, p_tot, p_1us, app = _load_pdf(pmt)

    dt_max = min(np.percentile(dt_us, 99.5), 5.5)
    area_max = min(np.percentile(area_pe, 99.5), 30.0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex="col")
    fig.suptitle(f"{pmt} — {label} — {len(dt_us)} afterpulses", fontsize=13)

    _panel(axes[0, 0], dt_us, area_pe, False,
           "归一化前 (raw counts)", dt_max, area_max)
    _dt_projection(axes[1, 0], dt_us, dt_max, "Δt 分布 (归一化前)")

    _panel(axes[0, 1], dt_us, area_pe, True,
           "归一化后 (PDF, ∫=1)", dt_max, area_max)
    _dt_projection(axes[1, 1], dt_us, dt_max, "Δt 分布 (归一化后)")

    axes[0, 1].text(0.02, 0.98,
                    f"APP={app:.4f}\np_ap_total={p_tot:.5f}\np_ap_1us={p_1us:.5f}",
                    transform=axes[0, 1].transAxes, va="top", fontsize=9,
                    bbox=dict(boxstyle="round", fc="w", alpha=0.7))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / f"app_pdf_{pmt}_2d_before_after.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pdf_marginals(pmt: str, label: str):
    """每 PMT 一张大图: 归一化 PDF 的时间边际 + 电荷边际 (两子画布)。"""
    h, dt_edges_us, npe_edges, p_tot, p_1us, app = _load_pdf(pmt)
    pdf = h / h.sum() if h.sum() else h

    dt_c = 0.5 * (dt_edges_us[:-1] + dt_edges_us[1:])
    npe_c = 0.5 * (npe_edges[:-1] + npe_edges[1:])

    # 边际分布 (归一化 PDF 在时间/电荷轴上的投影, 均为概率密度)
    pdf_dt = pdf.sum(axis=1)
    pdf_npe = pdf.sum(axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    fig.suptitle(f"{pmt} — {label} — 归一化 PDF 边际分布", fontsize=13)

    ax = axes[0]
    ax.plot(dt_c, pdf_dt, lw=1.2, color="steelblue")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Delta Time [$\mu$s]")
    ax.set_ylabel("Probability / bin")
    ax.set_title("时间边际分布 P(Δt)", fontsize=10)

    ax = axes[1]
    ax.step(npe_c, pdf_npe, where="mid", lw=1.2, color="seagreen")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Area [PE]")
    ax.set_ylabel("Probability / bin")
    ax.set_title("电荷边际分布 P(area)", fontsize=10)

    axes[0].text(0.03, 0.97, f"APP={app:.4f}\np_ap_total={p_tot:.5f}\np_ap_1us={p_1us:.5f}",
                 transform=axes[0].transAxes, va="top", fontsize=9,
                 bbox=dict(boxstyle="round", fc="w", alpha=0.7))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_DIR / f"app_pdf_{pmt}_marginals.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_delta_time_grid():
    """PMT Δt 分布网格, 时间单位 μs (风格同 app.py)。"""
    n_pmt = len(PMTS)
    cols = max(n_pmt, 1)
    rows = 1
    fig, axes = plt.subplots(rows, cols, figsize=(2.5 * cols, 2.2 * rows),
                             sharex=True, sharey=True)
    if n_pmt == 1:
        axes = np.array([[axes]])

    for idx, pmt in enumerate(PMTS):
        ax = axes[idx] if cols > 1 else axes[0, 0]
        dt_us, _ = _scatter(pmt)
        ax.hist(dt_us, bins=80, range=(0, 4),
                histtype="stepfilled", color="skyblue", edgecolor="skyblue",
                linewidth=0.8, alpha=0.9)
        ax.set_title(f"{pmt} ({PMTS[pmt][2]}) — {len(dt_us)} AP", fontsize=8)
        ax.set_xlim(0, 4)
        ax.tick_params(labelsize=7, length=3, width=0.8)
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci", axis="y", scilimits=(4, 4))
        ax.yaxis.get_offset_text().set_fontsize(6)

    for col in range(cols):
        ax = axes[col] if cols > 1 else axes[0, 0]
        ax.set_xlabel("Delay Time [$\mu$s]", fontsize=8)
    axes[0].set_ylabel("Counts", fontsize=8)

    fig.suptitle("Afterpulse Delta Time Distribution (typical PMTs)", fontsize=9, y=1.01)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(OUT_DIR / "app_pdf_delta_time_grid.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_2d_full_range(pmt: str, label: str, dt_max: float = 60.0, area_max: float = 40.0):
    """LV2358 原始散点二维图: dt 至 dt_max μs, area 至 area_max PE (log 色标)。"""
    dt_us, area_pe = _scatter(pmt)

    fig, ax = plt.subplots(figsize=(11, 7))
    hist = ax.hist2d(
        dt_us, area_pe,
        bins=[200, 100],
        range=[[0, dt_max], [0, area_max]],
        cmap="jet",
        density=True,
        norm=matplotlib.colors.LogNorm(),
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.12)
    cbar = plt.colorbar(hist[3], ax=cax)
    cbar.set_label(r"Density [arb.]", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlim(0, dt_max)
    ax.set_ylim(0, area_max)
    ax.set_xlabel(r"Time Delay [$\mu$s]", fontsize=12)
    ax.set_ylabel("Afterpulse Area [PE]", fontsize=12)
    ax.set_title(f"{pmt} — {label} — {len(dt_us)} afterpulses\n"
                 r"Afterpulse Area vs $\Delta$t (full range)", fontsize=12)
    ax.tick_params(direction="in", labelsize=9, length=4, width=1)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"app_{pmt}_2d_dt60_area40.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for pmt in PMTS:
        plot_2d_before_after(pmt, PMTS[pmt][2])
        plot_pdf_marginals(pmt, PMTS[pmt][2])
    plot_delta_time_grid()
    # 全范围二维图: 本次只出 LV2358 (dt→60us, area→40PE)
    plot_2d_full_range("LV2358", PMTS["LV2358"][2])
    print(f"plots written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
