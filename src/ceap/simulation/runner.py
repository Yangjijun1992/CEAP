"""模拟主流程 (集成 SIM 各模块 + 输出指标 S-05)。

流程: 加载配置 -> 构建 PDF/探测器/μ 子模型 -> 时间线生成(西蒙背景)
      -> μ 子 S1 主脉冲 + after pulse 叠加 -> 窗口扫描 -> 输出指标。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..config.loader import Config, load_config
from ..models.afterpulse_pdf import AfterpulsePDFFactory
from ..models.detector import Detector
from ..models.muon import MuonModel
from ..simulation.afterpulse_gen import AfterpulseGenerator
from ..simulation.timeline import TimelineSimulator
from ..simulation.window_scanner import WindowScanner


def run_simulation(cfg: Config):
    """执行一次完整模拟，返回结果字典并写入 output_dir。"""
    sim = cfg.get("simulation") or {}
    duration_us = float(sim.get("duration_us", 1e6))
    seed = sim.get("seed")
    out_dir = Path(sim.get("output_dir", "output"))
    run_id = sim.get("run_id", "run0")
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 构建模型 -------------------------------------------------
    pdf = AfterpulsePDFFactory.build(cfg)
    det = Detector(cfg)
    n_pmt = max(det.n_pmts, 1)
    muon = MuonModel(cfg, seed=seed)

    # --- 1. 时间线: 物理本底 + 暗噪声 ---------------------------
    tl = TimelineSimulator(cfg, seed=seed, n_pmt=n_pmt)
    hits, base_mains = tl.generate(duration_us)

    # --- 2. μ 子主脉冲: 到临时间 + S1 PE -----------------------
    rng = np.random.default_rng(seed)
    muon_mains = []
    if muon.hit_rate_hz > 0:
        n_mu = int(rng.poisson(muon.hit_rate_hz * duration_us * 1e-6))
        t_end_ns = duration_us * 1e3
        mt = rng.uniform(0, t_end_ns, size=n_mu)
        for t in mt:
            muon_mains.append(_MainPulseLike(t_ns=float(t), npe=muon.sample_main_s1_npe(), kind="muon"))
    main_pulses = list(base_mains) + list(muon_mains)
    main_pulses.sort(key=lambda m: m.t_ns)

    # --- 3. after pulse 叠加 (SIM-02) ----------------------------
    gen = AfterpulseGenerator(cfg, pdf, seed=seed, n_pmt=n_pmt)
    for m in main_pulses:
        hits.extend(gen.generate_for(m, n_pmt))

    # --- 4. 窗口扫描 (SIM-03) + 指标 (S-05) -------------------
    scanner = WindowScanner(cfg)
    events, n_windows, n_events = scanner.scan(hits, main_pulses, duration_us)
    rate_hz = scanner.compute_rate(n_events, duration_us)

    # --- 5. 汇总输出 --------------------------------------------
    result = {
        "run_id": run_id,
        "duration_us": duration_us,
        "n_pmt": n_pmt,
        "n_dark_physics_hits": len(hits),
        "n_main_pulses": len(main_pulses),
        "n_muon_mains": len(muon_mains),
        "muon_hit_rate_hz": muon.hit_rate_hz,
        "n_windows_scanned": n_windows,
        "n_trigger_events": n_events,
        "background_rate_hz": rate_hz,
        "rate_per_muon": (rate_hz / muon.hit_rate_hz) if muon.hit_rate_hz else 0.0,
        "pe_spectrum": _pe_spectrum(events),
        "pdf": pdf.describe(),
    }
    _write_outputs(out_dir, run_id, result, events, sim.get("save_windows", True))
    return result


def _MainPulseLike(t_ns, npe, kind):
    from ..simulation.timeline import MainPulse
    return MainPulse(t_ns=t_ns, npe=npe, kind=kind)


def _pe_spectrum(events):
    if not events:
        return []
    vals = np.array([e.total_pe for e in events])
    hist, edges = np.histogram(vals, bins=20)
    return {"counts": hist.tolist(), "edges": edges.tolist()}


def _write_outputs(out_dir: Path, run_id: str, result: dict, events, save_windows: bool):
    with open(out_dir / f"{run_id}_summary.json", "w") as f:
        json.dump(result, f, indent=2, default=str)
    if save_windows and events:
        arr = np.array([[e.t_ns, e.total_pe, int(e.has_physics_s1)] for e in events])
        np.save(out_dir / f"{run_id}_windows.npy", arr)
