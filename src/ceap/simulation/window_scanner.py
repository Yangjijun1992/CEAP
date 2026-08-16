"""1 μs 偶然符合窗口扫描 (SIM-03) 与输出指标 (S-05)。

在击中型列表上滑动/步进窗口，计算窗内总 PE，记录超阈值事例；
并统计噪声事例率、PE 谱。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..config.loader import Config
from ..simulation.timeline import Hit, MainPulse


@dataclass
class WindowEvent:
    t_ns: float
    total_pe: float
    has_physics_s1: bool


class WindowScanner:
    """对给定时长内的击中列表执行窗口扫描。"""

    def __init__(self, cfg: Config):
        sw = cfg.get("signal_window") or {}
        self.length_us = float(sw.get("length_us", 1.0))
        self.step_us = float(sw.get("step_us", 0.0))  # 0 表示连续滑动
        self.thr_min = float(sw.get("threshold_pe_min", 0.0))
        self.thr_max = float(sw.get("threshold_pe_max", np.inf))
        self.require_no_physics = bool(sw.get("require_no_physics_s1", True))
        dz = sw.get("dead_zone") or {}
        self.dead_s1_ns = float(dz.get("s1_after_ns", 0.0))

    def _is_in_dead_zone(self, t_ns: float, mains: list[MainPulse]) -> bool:
        """S-04: 检查某时刻是否落在某个主脉冲 S1 死区内。"""
        if not mains:
            return False
        for m in mains:
            if 0.0 <= t_ns - m.t_ns <= self.dead_s1_ns:
                return True
            if t_ns < m.t_ns:
                break
        return False

    def scan(self, hits: list[Hit], mains: list[MainPulse], duration_us: float):
        """执行滑动窗口扫描。

        返回 (events: list[WindowEvent], n_windows_long, n_events_long)。
        """
        mains_sorted = sorted(mains, key=lambda m: m.t_ns)
        times = np.array([h.t_ns for h in hits]) if hits else np.empty(0)
        npe = np.array([h.npe for h in hits]) if hits else np.empty(0)

        events: list[WindowEvent] = []
        length_ns = self.length_us * 1e3
        step_ns = self.step_us * 1e3 if self.step_us > 0 else 1.0
        t = 0.0
        t_end_ns = duration_us * 1e3
        n_long = 0
        n_events_long = 0

        while t + length_ns <= t_end_ns:
            lo, hi = t, t + length_ns
            inwin = (times >= lo) & (times < hi)
            total = float(npe[inwin].sum())
            has_phy = self._has_physics_s1(lo, hi, mains_sorted)
            triggered = self.thr_min <= total <= self.thr_max
            if not self.require_no_physics or not has_phy:
                if triggered:
                    events.append(WindowEvent(t_ns=lo + length_ns / 2, total_pe=total, has_physics_s1=has_phy))
            n_long += 1
            if triggered and not has_phy:
                n_events_long += 1
            t += step_ns
        return events, n_long, n_events_long

    @staticmethod
    def _has_physics_s1(lo_ns: float, hi_ns: float, mains: list[MainPulse]) -> bool:
        """窗口内是否包含物理 S1 主脉冲 (S-03)。"""
        for m in mains:
            if lo_ns <= m.t_ns < hi_ns:
                return True
            if m.t_ns >= hi_ns:
                break
        return False

    def compute_rate(self, n_events: int, duration_us: float) -> float:
        """估算噪声/本底事例率 (Hz)。"""
        return n_events / (duration_us * 1e-6) if duration_us > 0 else 0.0
