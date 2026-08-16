"""Afterpulse (时间差, 电荷) 二维 PDF 接口与实现。

核心需求 (AP-01..AP-06)：
  输入大信号后 after pulse 的 (时间差, 电荷大小) 二维 PDF。

数据来源：
  1) 外部文件 (pdf_file): npz 保存二维直方图 h[bin_dt, bin_npe]，以及
     时间轴 dt_edges[ns]、电荷轴 npe_edges、归一化参数
     p_ap_total（每主 PE 全程概率）、p_ap_1us（1 μs 窗口概率）、app。
     由 tools/app_pdf_builder.py 生成。
  2) 逐 PMT PDF (pmt_pdfs): {pmt_index: pdf_file}，见 PerPMTPDF (AP-05)。
  3) 参数化形式 (params): 占位，连接具体函数后替换实现。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from ..config.loader import Config


@dataclass
class BaseAfterpulsePDF:
    """Afterpulse PDF 的通用接口与公共方法。

    线性模型 (AP-02)：每个主脉冲 PE 独立以 p_ap_total 概率产生 after pulse。
    """

    seed: Optional[int] = None
    _rng: np.random.Generator = field(default=None, repr=False)

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)

    # -- 需子类实现 ---------------------------------------------
    def sample_once(self, rng: np.random.Generator):
        """采样 1 个 after pulse 的 (dt[ns], npe)。子类实现。"""
        raise NotImplementedError

    # -- 公共方法 -----------------------------------------------
    def sample(self, main_npe: int):
        """为一个主脉冲(main_npe 个 PE)生成整条 after pulse 序列。

        线性模型：每个 PE 独立以 p_ap_total 概率产生 after pulse。
        当 PDF 依赖主脉冲 PE 数时，子类应重写此方法。

        返回: np.ndarray shape (N, 2) 每行 (dt_ns, npe)
        """
        if main_npe <= 0 or self.p_ap_total <= 0:
            return np.empty((0, 2), dtype=float)
        n_events = self._rng.binomial(main_npe, min(self.p_ap_total, 1.0))
        if n_events == 0:
            return np.empty((0, 2), dtype=float)
        seq = np.empty((n_events, 2), dtype=float)
        for i in range(n_events):
            dt, npe = self.sample_once(self._rng)
            seq[i] = (dt, npe)
        return seq

    def p_in_window(self, dt_ns: float):
        """AP-03: dt 在 [0, dt_ns] 内的 after pulse 概率（每主 PE）。"""
        raise NotImplementedError

    def describe(self) -> str:
        return "BaseAfterpulsePDF (unconfigured)"


class FileAfterpulsePDF(BaseAfterpulsePDF):
    """从二维直方图文件加载的 after pulse PDF。

    npz 支持两套 key：
      新格式 (app_pdf_builder 输出): h, dt_edges, npe_edges, p_ap_total,
                                    p_ap_1us, app
      旧格式: h, dt_edges, npe_edges（p_ap_total 退化为 h 积分，不推荐）
    """

    def __init__(self, pdf_file: str, max_npe: int = 30, seed: Optional[int] = None, **kwargs):
        super().__init__(seed=seed)
        self.pdf_file = Path(pdf_file)
        self.max_npe = int(max_npe)
        self.load()

    def load(self):
        data = np.load(self.pdf_file)
        self.hist = np.asarray(data["h"], dtype=float)
        self.dt_edges = np.asarray(data["dt_edges"], dtype=float)   # ns
        self.npe_edges = np.asarray(data["npe_edges"], dtype=float)
        self.dt_centers = 0.5 * (self.dt_edges[:-1] + self.dt_edges[1:])
        self.npe_centers = 0.5 * (self.npe_edges[:-1] + self.npe_edges[1:])

        # 归一化: h 视为 (dt,npe) 联合事件计数
        self._n_events = float(self.hist.sum())
        self._pdf = self.hist / (self._n_events or 1.0)
        self._flat_pdf = self._pdf.ravel()
        self._flat_idx = np.flatnonzero(self._flat_pdf)

        # 归一化参数 (优先读文件, 旧格式退化)
        if "p_ap_total" in data.files:
            self.p_ap_total = float(data["p_ap_total"])
        else:
            self.p_ap_total = float(self._n_events)
        if "p_ap_1us" in data.files:
            self.p_ap_1us = float(data["p_ap_1us"])
        else:
            self.p_ap_1us = self.p_in_window(self._window_ns())
        if "app" in data.files:
            self.app = float(data["app"])
        else:
            self.app = float("nan")
        self.pmt_id = str(data.get("pmt_id", "unknown")) if "pmt_id" in data.files else "unknown"

    def _window_ns(self) -> float:
        return 1000.0

    # -- AP-04: max_npe 强制 ------------------------------------
    def sample_once(self, rng: np.random.Generator):
        idx = rng.choice(self._flat_idx, p=self._flat_pdf[self._flat_idx])
        i, j = np.unravel_index(idx, self.hist.shape)
        npe = float(self.npe_centers[j])
        # 强制电荷上限 (AP-04)
        npe = min(npe, float(self.max_npe))
        return float(self.dt_centers[i]), npe

    def p_in_window(self, dt_ns: float):
        """AP-03: 每主 PE 在 [0, dt_ns] 时间窗内产生 after pulse 的概率。

        按事件面积占比折算: 窗内事件面积 / 总面积 × p_ap_total。
        """
        dt_ns = float(dt_ns)
        mask = self.dt_centers < dt_ns
        if not mask.any():
            return 0.0
        area_in = float((self._pdf[mask, :] * self.npe_centers[None, :]).sum())
        area_all = float((self._pdf * self.npe_centers[None, :]).sum())
        if area_all <= 0:
            return 0.0
        return self.p_ap_total * (area_in / area_all)

    def describe(self) -> str:
        return (f"FileAfterpulsePDF({self.pdf_file.name}, bins={self.hist.shape}, "
                f"p_ap_total={self.p_ap_total:.4g}, p_ap_1us={self.p_ap_1us:.4g}, "
                f"app={self.app:.4g})")


class PerPMTPDF(BaseAfterpulsePDF):
    """逐 PMT PDF 集合 (AP-05)：每只 PMT 独立 PDF，支持分组平均。

    构造: {pmt_index: FileAfterpulsePDF}；未配置的 PMT 用默认 PDF。
    """

    def __init__(self, pdfs: dict[int, BaseAfterpulsePDF],
                 default_pdf: Optional[BaseAfterpulsePDF] = None, seed: Optional[int] = None):
        super().__init__(seed=seed)
        self.pdfs = {int(k): v for k, v in pdfs.items()}
        self.default_pdf = default_pdf

    def get(self, pmt_index: int) -> BaseAfterpulsePDF:
        return self.pdfs.get(int(pmt_index), self.default_pdf)

    def sample(self, main_npe: int):
        raise NotImplementedError("PerPMTPDF: 请使用 get(pmt).sample(npe) 逐 PMT 采样")

    def p_in_window(self, dt_ns: float):
        return 0.0

    def describe(self) -> str:
        parts = ", ".join(f"{k}:{v.pmt_id if hasattr(v, 'pmt_id') else '?'}" for k, v in sorted(self.pdfs.items()))
        return f"PerPMTPDF({{{parts}}}, n_pmt={len(self.pdfs)})"


class PlaceholderAfterpulsePDF(BaseAfterpulsePDF):
    """占位实现：默认不产生 after pulse，用于连接 PDF 前的流程调试。"""

    def __init__(self, **kwargs):
        super().__init__(seed=kwargs.get("seed"))
        self.max_npe = int(kwargs.get("max_npe", 1))
        self.p_ap_total = 0.0
        self.p_ap_1us = 0.0
        self.app = float("nan")
        self.pmt_id = "placeholder"

    def sample_once(self, rng: np.random.Generator):
        return (np.inf, 0.0)

    def p_in_window(self, dt_ns: float):
        return 0.0

    def describe(self) -> str:
        return "PlaceholderAfterpulsePDF (no afterpulse produced)"


class AfterpulsePDFFactory:
    """根据配置构建 after pulse PDF 对象。"""

    @staticmethod
    def build(cfg: Config) -> BaseAfterpulsePDF:
        ap = cfg.get("afterpulse") or {}
        seed = (cfg.get("simulation") or {}).get("seed")
        max_npe = ap.get("max_npe", 30)

        # AP-05: 逐 PMT PDF
        pmt_pdfs = ap.get("pmt_pdfs") or {}
        if pmt_pdfs:
            built = {}
            for pmt, f in pmt_pdfs.items():
                built[int(pmt)] = FileAfterpulsePDF(f, max_npe=max_npe, seed=seed)
            default = None
            if ap.get("pdf_file"):
                default = FileAfterpulsePDF(ap["pdf_file"], max_npe=max_npe, seed=seed)
            return PerPMTPDF(built, default_pdf=default, seed=seed)

        pdf_file = ap.get("pdf_file")
        if pdf_file:
            return FileAfterpulsePDF(pdf_file, max_npe=max_npe, seed=seed)
        return PlaceholderAfterpulsePDF(max_npe=max_npe, seed=seed)
