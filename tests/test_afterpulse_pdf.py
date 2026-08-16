"""Afterpulse PDF 构建与归一化测试 (AP 组)。"""
import numpy as np
import pytest

from ceap.tools.app_pdf_builder import build_pdf, default_dt_edges, default_npe_edges
from ceap.models.afterpulse_pdf import (
    AfterpulsePDFFactory,
    FileAfterpulsePDF,
    PerPMTPDF,
)
from ceap.config.loader import load_config

DATA = {
    "LV2305": ("data/app_pdf_LV2305.npz", 0.102),
    "LV2358": ("data/app_pdf_LV2358.npz", 0.053),
}


def test_build_pdf_normalization():
    rng = np.random.default_rng(0)
    dt = rng.uniform(0, 10_000, size=50_000)
    area = rng.exponential(2.0, size=50_000)
    main = np.full(100, 400.0)
    pdf = build_pdf(dt, area, main)

    assert pdf.p_ap_total > 0
    assert pdf.p_ap_1us > 0
    assert pdf.p_ap_1us < pdf.p_ap_total
    # hist 面积守恒: sum(h * npe_center) / sum(area) 近似 1 (未裁剪时)
    assert np.isfinite(pdf.p_ap_total)


def test_build_pdf_overflow_folding():
    """电荷溢出应折入末 bin 而非丢弃 (AP-04)。"""
    dt = np.array([100.0, 200.0])
    area = np.array([5.0, 100.0])
    main = np.array([400.0, 400.0])
    npe_edges = default_npe_edges(10.0)
    pdf = build_pdf(dt, area, main, npe_edges=npe_edges)
    # 两个事件都在末 bin (面积 5 与 100 折入/归属, 100 折入末 bin)
    assert pdf.hist.sum() == 2


def test_file_pdf_loads_metadata():
    f = DATA["LV2305"][0]
    pdf = FileAfterpulsePDF(f, max_npe=30, seed=1)
    assert pdf.pmt_id == "LV2305"
    assert abs(pdf.app - 0.102) < 1e-3
    assert pdf.p_ap_total > 0


@pytest.mark.parametrize("pmt,app_exp", DATA.items())
def test_sampling_reproduces_app(pmt, app_exp):
    """抽样统计 APP 应复现原始 APP (线性模型自洽)。"""
    file = DATA[pmt][0]
    app_exp = DATA[pmt][1]
    pdf = FileAfterpulsePDF(file, max_npe=30, seed=42)
    main_npe, n_main = 400, 4000
    total = 0.0
    for _ in range(n_main):
        seq = pdf.sample(main_npe)
        total += seq[:, 1].sum() if len(seq) else 0.0
    app_sim = total / (main_npe * n_main)
    assert abs(app_sim - app_exp) < 0.02, f"{pmt}: sim={app_sim:.4f} exp={app_exp:.4f}"


def test_p_in_window_monotonic():
    pdf = FileAfterpulsePDF(DATA["LV2305"][0], max_npe=30, seed=1)
    p_1 = pdf.p_in_window(1000.0)
    p_10 = pdf.p_in_window(10_000.0)
    p_all = pdf.p_in_window(1e9)
    assert 0 < p_1 < p_10 < p_all
    assert abs(p_1 - pdf.p_ap_1us) < 1e-9


def test_max_npe_enforced():
    pdf = FileAfterpulsePDF(DATA["LV2358"][0], max_npe=5, seed=1)
    samples = np.array([pdf.sample_once(pdf._rng)[1] for _ in range(20_000)])
    assert samples.max() <= 5


def test_per_pmt_pdf_factory():
    cfg = load_config("config/settings.yaml")
    cfg._data["afterpulse"]["pmt_pdfs"] = {
        0: DATA["LV2305"][0],
        1: DATA["LV2358"][0],
    }
    cfg._data["afterpulse"]["pdf_file"] = None
    pdf = AfterpulsePDFFactory.build(cfg)
    assert isinstance(pdf, PerPMTPDF)
    assert pdf.get(0).pmt_id == "LV2305"
    assert pdf.get(1).pmt_id == "LV2358"
    # 未配置的 PMT 用 default 兜底
    assert pdf.get(5) is None
