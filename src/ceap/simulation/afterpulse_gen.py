"""Afterpulse 生成器 (SIM-02)。

对每个主脉冲(含 μ 子 S1 与物理本底)按 PDF 独立为各 PMT 生成
after pulse 的时间与电荷。

支持两种 PDF 形态:
  - 单一 PDF: 所有 PMT 共用 (uniform, AP-05)
  - PerPMTPDF: 每只 PMT 独立 PDF (per_pmt 分组, AP-05)
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..config.loader import Config
from ..models.afterpulse_pdf import BaseAfterpulsePDF, PerPMTPDF
from ..simulation.timeline import Hit, MainPulse


class AfterpulseGenerator:
    """根据主脉冲 PE 数与 PDF 产生 after pulse 击中型。"""

    def __init__(self, cfg: Config, pdf: BaseAfterpulsePDF, seed: Optional[int] = None,
                 n_pmt: int = 1):
        self.cfg = cfg
        self.pdf = pdf
        self.n_pmt = n_pmt
        self._rng = np.random.default_rng(seed)

    def generate_for(self, main: MainPulse, n_pmt: int) -> list[Hit]:
        """为一个主脉冲生成全 PMT 的 after pulse 序列。

        简化线性模型：总 npe 平均分配到各 PMT(占位，T-M-04 细化后替换)，
        每只 PMT 使用其对应的 PDF 独立抽样。

        返回 Hit 列表。
        """
        out: list[Hit] = []
        if main.npe <= 0:
            return out
        per_pmt = max(int(main.npe // max(n_pmt, 1)), 1)
        for pmt in range(n_pmt):
            pdf = self.pdf.get(pmt) if isinstance(self.pdf, PerPMTPDF) else self.pdf
            seq = pdf.sample(per_pmt)
            for dt_ns, npe in seq:
                if np.isfinite(dt_ns):
                    out.append(Hit(pmt=pmt, t_ns=main.t_ns + dt_ns, npe=npe))
        return out
