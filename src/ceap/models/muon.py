"""Muon model (M 组): 生成 μ 子 S1 主脉冲的 (时间, PE 数)。

需求:
  M-01 通量/角分布  M-02 几何接收度  M-03 能损模型
  M-04 S1 产额与分配 M-05 校验样本(TBD)
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..config.loader import Config


class MuonModel:
    """μ 子模型：返回 μ 子到临时间与其主脉冲总 PE 数。"""

    def __init__(self, cfg: Config, seed: Optional[int] = None):
        self._rng = np.random.default_rng(seed)
        mu = cfg.get("muon") or {}
        self.flux_per_m2_hz = float(mu.get("flux_per_m2_hz", 0.0))
        self.active_area_m2 = float(mu.get("active_area_m2", 0.0))
        self.dedx_mev_cm = float(mu.get("dedx_mev_cm", 2.0))
        track_mean_cm = float(mu.get("track_length_mean_cm", 0.0))
        self.track_length_mean_cm = track_mean_cm
        # 光产额：优先用 muon 里的覆盖值，否则用 detector 值
        ly = mu.get("light_yield_pe_keV")
        self.ly = float(ly if ly else (cfg.get("detector") or {}).get("target", {}).get("light_yield_pe_keV", 0.0))
        self.s1_sigma_rel = float(mu.get("s1_sigma_rel", 0.1))
        # S1 时间区间 (用于死区): 取 detector.daq 波形长度或默认
        dag = (cfg.get("detector") or {}).get("daq") or {}
        self.s1_duration_ns = float(dag.get("waveform_length_us", 0.0)) * 1e3 or 1000.0
        self._hit_rate = self._compute_hit_rate()

    def _compute_hit_rate(self) -> float:
        """M-02: 几何接收度 × 通量 = 击中率 (Hz)。"""
        return self.flux_per_m2_hz * self.active_area_m2

    @property
    def hit_rate_hz(self) -> float:
        return self._hit_rate

    def mean_energy_mev(self) -> float:
        """M-03: 平均能量沉积 = dE/dx × 平均径迹长度。"""
        return self.dedx_mev_cm * self.track_length_mean_cm

    def sample_main_s1_npe(self) -> int:
        """M-04: 采样一个 μ 子主脉冲的总 PE 数。

        能量沉积 ~ 平均能损; S1 光子 -> PE, 加光子统计涨落 (简化)。
        """
        e_kev = self.mean_energy_mev() * 1e3
        npe_mean = e_kev * self.ly
        if npe_mean <= 0:
            return 0
        sigma = max(npe_mean * self.s1_sigma_rel, 1.0)
        return int(round(max(self._rng.normal(npe_mean, sigma), 0.0)))
