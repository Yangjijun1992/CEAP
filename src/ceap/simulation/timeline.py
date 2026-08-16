"""时间线模拟器 (SIM-01)。

按事例率向每个 PMT 注入物理背景 + 暗噪声，生成击中型列表，
并为 μ 子 S1 记录“主脉冲”(用于后续 after pulse 生成)。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..config.loader import Config


@dataclass
class Hit:
    pmt: int
    t_ns: float
    npe: float


@dataclass
class MainPulse:
    t_ns: float
    npe: int
    kind: str  # "muon" | "background"


class TimelineSimulator:
    """生成指定时长内的物理本底/暗噪声击中列表与主脉冲列表。"""

    def __init__(self, cfg: Config, seed: Optional[int] = None, n_pmt: int = 0):
        self.cfg = cfg
        self._rng = np.random.default_rng(seed)
        det = cfg.get("detector") or {}
        self.n_pmt = n_pmt or len((det.get("pmt_list") or [])) or 1
        # 单 PE 暗计数率 (每个 PMT, 单位 Hz) [TBD D-05]，可在 config 中指定
        tl = cfg.get("simulation", {}).get("timeline", {}) if isinstance(cfg.get("simulation"), dict) else {}
        self.dark_rate_hz = float(tl.get("dark_rate_hz_per_pmt", 0.0))
        # 物理本底 S1 触发率 (单位 Hz) [TBD D-05]
        self.background_rate_hz = float(tl.get("background_rate_hz", 0.0))

    def add_dark_noise(self, hits, t_start_ns, t_end_ns):
        """注入单 PE 暗噪声。返回新增 Hit 列表。"""
        if self.dark_rate_hz <= 0:
            return
        dt = (t_end_ns - t_start_ns) * 1e-9  # s
        n_total = self._rng.poisson(self.dark_rate_hz * self.n_pmt * dt)
        pmts = self._rng.integers(0, self.n_pmt, size=n_total)
        times = self._rng.uniform(t_start_ns, t_end_ns, size=n_total)
        for pmt, t in zip(pmts, times):
            hits.append(Hit(pmt=int(pmt), t_ns=float(t), npe=1.0))

    def add_background_mainpulses(self, mains, t_start_ns, t_end_ns):
        """注入物理本底主脉冲 (S-03: 计入无物理 S1 判断)。"""
        if self.background_rate_hz <= 0:
            return
        dt = (t_end_ns - t_start_ns) * 1e-9
        n = self._rng.poisson(self.background_rate_hz * dt)
        times = self._rng.uniform(t_start_ns, t_end_ns, size=n)
        for t in times:
            mains.append(MainPulse(t_ns=float(t), npe=int(self._rng.poisson(50)), kind="background"))

    def generate(self, duration_us: float):
        """生成 [0, duration_us] 内的击中列表与主脉冲列表。

        返回 (hits, mains)。
        """
        hits: list[Hit] = []
        mains: list[MainPulse] = []
        t_end_ns = duration_us * 1e3
        self.add_dark_noise(hits, 0.0, t_end_ns)
        self.add_background_mainpulses(mains, 0.0, t_end_ns)
        return hits, mains
