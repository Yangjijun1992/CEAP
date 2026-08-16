"""Detector descriptor (D 组参数容器)。

仅承载探测器参数，供模拟/分析模块读取。索引约定：
  - pmt_index: 0..n_pmts-1
"""
from __future__ import annotations

from ..config.loader import Config


class Detector:
    """探测器参数封装，从配置加载 [TBD] 占位字段。"""

    def __init__(self, cfg: Config):
        det = cfg.get("detector") or {}
        self.n_pmts = int(det.get("n_pmts", 0)) or len(det.get("pmt_list", []))
        self.pmt_list = list(det.get("pmt_list", []))
        self.single_pe = det.get("single_pe", {})
        self.target = det.get("target", {})
        self.daq = det.get("daq", {})
        self.background_rates = det.get("background_rates", [])
        self.environment = det.get("environment", {})

    @property
    def s1_duration_ns(self) -> float:
        return float(self.daq.get("waveform_length_us", 0.0)) * 1e3 or 1000.0
