# Two-Stream Instability Observatory
### NASA MMS Live Data Monitor + Interactive Theory Laboratory

**Vivek Shrivastav · Sikkim University · 2025**

[![MMS TSI Analyzer](https://github.com/Vivek-Shrivastav/Two_Stream_Lab/actions/workflows/analyze.yml/badge.svg)](https://github.com/Vivek-Shrivastav/Two_Stream_Lab/actions)
[![Live Dashboard](https://img.shields.io/badge/Dashboard-Live-brightgreen)](https://vivek-shrivastav.github.io/Two_Stream_Lab/)

---

## What This Project Does

This is a **live research tool** that combines two things:

1. **Interactive Theory Laboratory** — the full mathematical framework for two-stream instability (TSI), including exact dispersion equations (Eq.18–21), growth rate calculators, SMILEI PIC simulation data, phase-space animation, and an experiment playground.

2. **Live MMS Data Observatory** — downloads real NASA MMS-1 spacecraft data every 6 hours and tests whether the theoretical conditions for TSI are actually met in the magnetosphere. Detects electron holes, beam activity, temperature anisotropy, and Langmuir waves. Maintains a 6-month rolling history and logs unusual plasma events.

---

## The Physics We Study

### Two-Stream Instability (TSI)
When two electron beams travel in opposite directions through each other (or one beam through a stationary background), the configuration is unstable if:

```
∂f/∂v > 0  at  v_phase = ω_pe / k
```

This positive slope in the velocity distribution function (bump-on-tail) allows the wave to extract energy from the beam, growing exponentially.

**Growth rates from our paper (Shrivastav et al. 2025 (https://arxiv.org/abs/2508.14362)):**

| Equation | Regime | Formula |
|---|---|---|
| Eq.(18) | Non-relativistic, single-species | γ = √{ ½[2v₀²k² + ω²ₚ − √(8v₀²k²ω²ₚ + ω⁴ₚ)] } |
| Eq.(19) | Relativistic, single-species | γ = √{ ½[2v₀²k² + ω²ₚ/γ₀³ − √(8v₀²k²ω²ₚ/γ₀³ + ω⁴ₚ/γ₀⁶)] } |
| Eq.(20) | Non-relativistic, multi-species | γ = √{ ½[2v₀²k² + 2ω²ₚ − √(16v₀²k²ω²ₚ + 4ω⁴ₚ)] } |
| Eq.(21) | Relativistic, multi-species | γ = √{ ½[2v₀²k² + 2ω²ₚ/γ₀³ − √(16v₀²k²ω²ₚ/γ₀³ + 4ω⁴ₚ/γ₀⁶)] } |

where γ₀ = 1/√(1−v₀²) is the Lorentz factor, v₀ is the beam velocity (in units of c), and ωₚ = 1 in normalised units.

All equations are validated against SMILEI PIC simulations (data embedded in the page).

### What We Look For in MMS Data

| Signature | What It Means | Detection Method |
|---|---|---|
| Electron holes | Nonlinear TSI saturation products — beam trapped particles | Bipolar E∥ pulses in EDP |
| Beam (df/dv > 0) | Classical bump-on-tail condition active | V_e >> v_thermal (FPI-DES) |
| TSI→Weibel chain | TSI scatters beam → builds T⊥ > T∥ → Weibel | Cross-correlation V_e vs T⊥/T∥ |
| Langmuir waves | Linear phase wave growth at f_pe | E-field PSD enhancement near f_pe |

### TSI → Weibel Energy Chain
PIC simulations (Innocenti et al. 2017) predict:
```
Electron beam  →  Langmuir waves grow  →  beam electrons scattered ⊥
→  T_perp > T_par builds  →  Weibel instability drives B fluctuations
→  Electromagnetic turbulence cascade
```
This causal chain has **not yet been directly confirmed observationally**. The cross-correlation analysis in this tool tests for it in every 6-hour MMS window.

---

## Repository Structure

```
Two_Stream_Instability/
│
├── index.html              ← Complete dashboard (theory + observatory)
├── fetcher.py              ← MMS data downloader and analyser
│
├── results/
│   ├── mms_live.json       ← Latest 6-hour window analysis
│   ├── mms_history.json    ← 6-month rolling history (max 720 records)
│   └── mms_events.json     ← Unusual event log (max 500 events)
│
└── .github/workflows/
    └── analyze.yml         ← GitHub Actions scheduler (every 6 hours)
```

---

## How It Works

### Data Pipeline

```
GitHub Actions triggers (every 6 hours)
         ↓
fetcher.py downloads MMS L2 data via NASA HAPI API
         ↓
         ├─ FGM: Magnetic field B (8 sps survey) ─────────────── Turbulence PSD, spectral index α
         ├─ EDP: Electric field E (32 sps fast) ──────────────── Electron hole detection, Langmuir PSD
         ├─ FPI-DES: Electron moments (~4.5 s) ────────────────── n_e, V_e, T_par, T_perp, f_pe
         └─ FPI-DIS: Ion moments (~4.5 s) ─────────────────────── n_i (for β, v_A)
         ↓
Derived analyses:
  ├─ Theory validation (Eq.18–21 with measured v₀, ω_pe)
  ├─ Electron hole detection (bipolar E∥ algorithm)
  ├─ df/dv slope proxy (v_beam/v_thermal)
  ├─ TSI→Weibel cross-correlation
  ├─ Langmuir wave proxy (PSD enhancement at f_pe)
  └─ Unusual event detection (6 categories)
         ↓
Results saved to results/*.json
         ↓
GitHub commits updated JSON files
         ↓
index.html reads JSON and renders live dashboard
```

### MMS Data Details

| Instrument | Dataset ID (HAPI) | Quantity | Rate |
|---|---|---|---|
| FGM | `MMS1_FGM_SRVY_L2@0` | B vector (GSE), nT | 8 sps |
| EDP | `MMS1_EDP_FAST_L2_DCE` | E vector (GSE), mV/m | 32 sps |
| FPI-DES | `MMS1_FPI_FAST_L2_DES-MOMS` | n_e, V_e, T_par, T_perp | ~4.5 s |
| FPI-DIS | `MMS1_FPI_FAST_L2_DIS-MOMS` | n_i, T_i | ~4.5 s |

**Important:** MMS L2 data has a ~90 day latency on NASA public servers. The fetcher automatically targets a window from 90 days ago to guarantee data availability.

---

## Unusual Event Categories

The fetcher automatically detects and logs 6 types of unusual events:

| Event Type | Trigger Condition | Physics Explanation |
|---|---|---|
| `HIGH_B_FIELD` | B_max > B_mean + 3σ, or B > 500 nT | Magnetospheric compression, current sheet, reconnection inflow |
| `ELECTRON_HOLE_BURST` | ≥ 10 bipolar E∥ events | Active/recent TSI nonlinear saturation phase |
| `HIGH_TEMPERATURE_ANISOTROPY` | A = T⊥/T∥ − 1 > 0.3 | Strong Weibel drive; beam-driven perpendicular heating |
| `ELECTRON_BEAM_DETECTED` | V_e > 300 km/s | TSI driver present; bump-on-tail condition may be met |
| `TSI_WEIBEL_CHAIN` | Beam + A > 0.05 simultaneously | Full TSI→Weibel energy cascade conditions present |
| `STEEP_TURBULENCE_SPECTRUM` | α < −2.5 | Enhanced kinetic-scale dissipation; Landau damping |

---

## TSI Activity Score

Each 6-hour window receives a score from 0–13:

| Evidence | Score |
|---|---|
| Electron beam detected (V_e > 300 km/s) | +2 |
| > 5 electron holes detected | +3 |
| Any electron holes detected | +1 |
| Langmuir wave enhancement | +2 |
| TSI→Weibel chain detected | +3 |
| Temperature anisotropy A > 0.05 | +1 |
| v_beam/v_thermal > 3 (bump-on-tail) | +2 |

| Score | Status |
|---|---|
| ≥ 7 | **ACTIVE** |
| 4–6 | **PROBABLE** |
| 2–3 | **POSSIBLE** |
| 0–1 | **QUIET** |

---

## Setup and Deployment

### Prerequisites
- GitHub account
- GitHub Pages enabled on your repository
- No local software needed (everything runs in GitHub Actions)

### Step 1: Create Repository
```bash
# Create a new repo: Two_Stream_Instability
# Upload: index.html, fetcher.py, results/*.json
```

### Step 2: Enable GitHub Pages
`Settings → Pages → Source: Deploy from branch → main → / (root)`

### Step 3: Set Up Workflow
Create `.github/workflows/analyze.yml`:
```yaml
name: MMS TSI Analyzer
on:
  schedule:
    - cron: '0 */6 * * *'
  workflow_dispatch:
permissions:
  contents: write
jobs:
  analyze:
    runs-on: ubuntu-latest
    timeout-minutes: 30
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: pip install numpy requests
      - name: Run MMS TSI analyzer
        run: python fetcher.py
      - name: Commit results
        run: |
          git config user.name  "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          git add results/
          git diff --staged --quiet || git commit -m "MMS TSI: $(TZ='Asia/Kolkata' date '+%d %b %Y %I:%M %p IST')"
          git push
```

### Step 4: Run First Fetch
`Actions → MMS TSI Analyzer → Run workflow`

The dashboard will populate after the first successful run (~30–60 seconds).

---

## Theory Laboratory Sections

| Section | Content |
|---|---|
| 01 | Physics concepts — beginner and expert mode |
| 02 | Growth rate calculator — exact Eq.18–21, single point and plot |
| 03 | Live growth rate lab — sliders, real-time γ(k) curves |
| 04 | Multi-curve comparison — all 4 regimes + SMILEI PIC data |
| 05 | Simulation box designer — SMILEI parameter calculator |
| 06 | Phase space explorer — live particle trapping animation |
| 07 | Experiment playground — save and overlay multiple runs |

---

## Observatory Sections

| Section | Content |
|---|---|
| 08 | Live MMS data — all 4 instruments, parameter table, B field plots |
| 09 | Theory validation — Eq.18–21 predictions vs measured conditions |
| 10 | Electron holes — bipolar E∥ detection, event table |
| 11 | Beam & df/dv — bump-on-tail test, V_e time series |
| 12 | TSI→Weibel chain — temperature anisotropy, cross-correlation |
| 13 | Langmuir waves — E-field PSD near f_pe |
| 14 | 6-month history — all quantities over time |
| 15 | Unusual events log — auto-detected anomalies with explanations |

---

## Key References

| Paper | Relevance |
|---|---|
| Shrivastav et al. (2025), Sikkim University | Equations 18–21, SMILEI PIC data in this tool |
| Buneman (1958) Phys. Rev. 115, 503 | Original two-stream instability derivation |
| Weibel (1959) Phys. Rev. Lett. 2, 83 | Weibel instability from temperature anisotropy |
| Innocenti et al. (2017) ApJL | TSI→Weibel energy chain (PIC simulations) |
| Bernstein, Greene, Kruskal (1957) PRL | Electron holes / BGK modes (nonlinear saturation) |
| Graham et al. (2016) GRL 43, 4098 | MMS electron hole observations |
| Steinvall et al. (2021) GRL | MMS bipolar E∥ structures |
| Jebaraj et al. (2025) A&A | Bump-on-tail vs truncated distribution debate |
| Malaspina et al. (2014) GRL 41, 5204 | MMS Langmuir wave observations |
| Burch et al. (2016) Science 352, aaf2939 | MMS mission description |
| Kolmogorov (1941) | Turbulence spectral index −5/3 |
| Iroshnikov (1964), Kraichnan (1965) | Turbulence spectral index −3/2 |

---

## Dependencies

```
Python packages:   numpy, requests   (pip install numpy requests)
JavaScript:        Plotly 2.27 (CDN), Google Fonts (CDN)
Data source:       NASA CDAS HAPI API (https://cdaweb.gsfc.nasa.gov/hapi)
```

No API keys required. All data is publicly available.

---

## Notes

- **Data latency:** MMS L2 survey data is published ~90 days after collection. The fetcher automatically uses the most recent available window (90 days ago).
- **f_pe and Langmuir detection:** The plasma frequency in the outer magnetosphere (f_pe ~ 1–100 kHz) is far above the EDP slow-survey Nyquist frequency (~16 Hz). Langmuir wave detection via PSD works best for high-density plasmas. For full Langmuir detection, EDP burst-mode data (8192 sps) would be needed.
- **History building:** The 6-month history builds up over time. After the first run you get 1 record; after 1 month ~120 records; after 6 months ~720 records.
- **SMILEI PIC validation:** The simulation data embedded in sections 04–05 covers v₀ = 0.1 to 0.95c. Agreement between analytical and PIC growth rates is typically within 5%.

---

*For questions or collaboration: Sikkim University, Department of Physics*
