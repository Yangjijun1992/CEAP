"""Afterpulse 二维 PDF 构建工具 (AP 组)。

将 `pmt_analysis` 仓库导出的原始散点数据
(ap_delta_time_ns, ap_area_pe) 转换为模拟框架使用的二维直方图 PDF，
并计算归一化参数:
  - p_ap_total: 主脉冲每 PE 在全部时间内的总 after pulse 概率 (AP-03)
  - p_ap_1us:   1 μs 窗口内概率 (AP-03)
  - APP:        原始定义 sum(ap_area_pe)/sum(main_area_pe)

归一化模型 (AP-02 线性模型):
  每个主脉冲 PE 独立产生 after pulse 事件, 事件概率 per-PE = APP / <area>.
  其中 <area> 为 after pulse 事件的面积(PE)期望。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class AppPdf:
    """构建后的二维 PDF 数据。

    Attributes:
        hist:  二维直方图, shape (n_dt_bins, n_npe_bins), 单位面积数(PE)
        dt_edges:  时间差轴边界 (ns)
        npe_edges: 电荷轴边界 (PE)
        p_ap_total: 每主 PE 全程 after pulse 概率
        p_ap_1us:   每主 PE 在 1 μs 窗口内概率
        app:        原始 APP 值
    """

    hist: np.ndarray
    dt_edges: np.ndarray
    npe_edges: np.ndarray
    p_ap_total: float
    p_ap_1us: float
    app: float

    def save(self, path: str | Path):
        np.savez(
            path,
            h=self.hist,
            dt_edges=self.dt_edges,
            npe_edges=self.npe_edges,
            p_ap_total=self.p_ap_total,
            p_ap_1us=self.p_ap_1us,
            app=self.app,
            pmt_id=getattr(self, "pmt_id", "unknown"),
        )


def default_dt_edges(max_dt_ns: float = 60_000.0, bin_ns: float = 100.0) -> np.ndarray:
    """默认时间轴: 0-60 μs, 步长 100 ns (覆盖离子峰区与长尾)。"""
    n = int(np.ceil(max_dt_ns / bin_ns))
    return np.linspace(0.0, n * bin_ns, n + 1)


def default_npe_edges(max_npe: float = 30.0) -> np.ndarray:
    """默认电荷轴: 0.5 起步, 步长 1 PE, 上限 max_npe (AP-04)。"""
    n = int(np.ceil(max_npe))
    return np.linspace(0.5, n + 0.5, n + 1)


def build_pdf(
    ap_delta_time_ns: np.ndarray,
    ap_area_pe: np.ndarray,
    main_area_pe: np.ndarray,
    dt_edges: np.ndarray | None = None,
    npe_edges: np.ndarray | None = None,
    window_us: float = 1.0,
) -> AppPdf:
    """从散点数据构建二维 PDF。

    Args:
        ap_delta_time_ns: after pulse 相对主脉冲结束的时间差 (ns)
        ap_area_pe: after pulse 面积 (PE)
        main_area_pe: 主脉冲面积 (PE)
        dt_edges: 时间轴边界 (ns), None 用默认
        npe_edges: 电荷轴边界 (PE), None 用默认
        window_us: p_ap_1us 的窗口长度 (μs)

    Returns:
        AppPdf: 归一化后的 PDF 及 p_ap_total / p_ap_1us。
    """
    dt = np.asarray(ap_delta_time_ns, dtype=float)
    area = np.asarray(ap_area_pe, dtype=float)
    main = np.asarray(main_area_pe, dtype=float)

    if dt_edges is None:
        dt_edges = default_dt_edges()
    if npe_edges is None:
        npe_edges = default_npe_edges()

    # 剔除无效值
    mask = np.isfinite(dt) & np.isfinite(area) & (dt >= 0) & (area > 0)
    dt, area = dt[mask], area[mask]

    # 溢出折入: 面积 > 最末 bin 的事件折入末 bin, 守恒 APP (AP-04)
    npe_max = float(npe_edges[-1])
    area_binned = np.clip(area, 0.0, npe_max - 1e-9)

    # 二维直方图: 面积数(PE), 非归一化
    hist, _, _ = np.histogram2d(dt, area_binned, bins=[dt_edges, npe_edges])
    n_events = int(hist.sum())

    # 事件面积期望 <area> (用 bin 中心, 与抽样一致)
    npe_centers = 0.5 * (npe_edges[:-1] + npe_edges[1:])
    area_mean = float((hist * npe_centers[None, :]).sum() / n_events) if n_events else 0.0

    # APP (原始定义) 与 per-PE 概率
    app = float(area.sum() / main.sum()) if main.size and main.sum() > 0 else 0.0
    # 线性模型: 每主 PE 产生 after pulse 事件的概率 = APP / <area>
    p_ap_total = app / area_mean if area_mean > 0 else 0.0

    # 1 μs 窗口内概率: 时间窗内事件面积占比 × per-PE 概率
    w_ns = window_us * 1e3
    dt_centers = 0.5 * (dt_edges[:-1] + dt_edges[1:])
    in_win = dt_centers < w_ns
    area_in_win = float((hist[in_win, :] * npe_centers[None, :]).sum())
    p_ap_1us = p_ap_total * (area_in_win / (area_mean * n_events)) if n_events else 0.0

    return AppPdf(hist=hist, dt_edges=dt_edges, npe_edges=npe_edges,
                  p_ap_total=p_ap_total, p_ap_1us=p_ap_1us, app=app)


def load_app_npz(path: str | Path) -> dict:
    """加载 pmt_analysis 导出的原始散点 npz。"""
    return dict(np.load(path, allow_pickle=True))


def build_from_npz(
    src: str | Path,
    dt_edges: np.ndarray | None = None,
    npe_edges: np.ndarray | None = None,
    window_us: float = 1.0,
    pmt_id: str | None = None,
) -> AppPdf:
    """从原始 npz 文件直接构建 PDF。"""
    data = load_app_npz(src)
    pdf = build_pdf(
        data["ap_delta_time_ns"],
        data["ap_area_pe"],
        data["main_area_pe"],
        dt_edges=dt_edges,
        npe_edges=npe_edges,
        window_us=window_us,
    )
    pdf.pmt_id = pmt_id or str(data.get("pmt_id", "unknown"))
    return pdf
