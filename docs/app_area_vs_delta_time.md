# Afterpulse Area vs Delta-time Distributions (Typical PMTs)

Selected typical PMT afterpulse distributions saved for later analysis. The
data come from the `pmt_analysis` repo's after-pulse (APP) analysis; each entry
is the per-event afterpulse **area (PE)** vs **delta time (ns after main
pulse)** scatter that feeds the APP 2D-histogram plots.

## Source repo

- Repo: `https://github.com/Yangjijun1992/pmt_analysis`
- APP analysis code: `src/pmt_analysis/analysis/app.py` (`analyze_app`,
  `find_afterpulse_candidates_per_channel`, `compute_app_per_channel`)
- Export script: `scripts/save_app_distributions.py`

## Selected PMTs (APP ≈ 20% / 10% / 5%)

| PMT | Run | APP (PE) | afterpulse count | main pulse count | data (.npz) | plot (.png) |
|-----|-----|----------|------------------|------------------|-------------|-------------|
| LV2229 | 00306 | 0.191 (~20%) | 74,700 | 1,897 | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00306_LV2229.npz` | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00306_LV2229.png` |
| LV2305 | 00381 | 0.102 (~10%) | 1,861,473 | 100,000 | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00381_LV2305.npz` | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00381_LV2305.png` |
| LV2358 | 00376 | 0.053 (~5%) | 1,197,530 | 100,000 | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00376_LV2358.npz` | `/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00376_LV2358.png` |

## Data directory (all files)

```
/mnt/data/PMT/R8520_406/output/app_distributions/
├── app_dist_00306_LV2229.npz
├── app_dist_00306_LV2229.png
├── app_dist_00381_LV2305.npz
├── app_dist_00381_LV2305.png
├── app_dist_00376_LV2358.npz
└── app_dist_00376_LV2358.png
```

## .npz contents (per file)

| key | description |
|-----|-------------|
| `run_id` | run id (str) |
| `pmt_id` | PMT id (str) |
| `channel` | channel index |
| `ap_delta_time_ns` | afterpulse delay vs main pulse end (ns), array |
| `ap_area_raw` | raw afterpulse area (ADC-integral scaled), array |
| `ap_area_pe` | PE-normalized afterpulse area, array |
| `main_area_pe` | main pulse area in PE, array |
| `app_pe` | channel APP = sum(ap_area_pe)/sum(main_area_pe) |

## Definition

```
APP = sum(afterpulse_area_pe) / sum(main_area_pe)
```

- Afterpulse search window: after `main_pulse.end + 35` samples, threshold
  `<-20 ADC`, min 2-sample dedup, then `findpulse_st_ed` per candidate.
- PE normalization uses per-PMT SPE gain (from the DB); for runs without
  gains the raw area is divided by the mean main-pulse area (fallback in the
  export script).
- Ion peaks in the delta-time axis: H⁺ 0.28 µs, He⁺ 0.56 µs, CH₄⁺ 1.01 µs,
  N₂⁺ 1.33 µs, Ar⁺ 1.58 µs, Xe⁺⁺ 2.02 µs, Xe⁺ 2.85 µs.

## Usage (read for later analysis)

```python
import numpy as np
d = np.load("/mnt/data/PMT/R8520_406/output/app_distributions/app_dist_00381_LV2305.npz")
dt_ns  = d["ap_delta_time_ns"]
area_pe = d["ap_area_pe"]
```
