# Contributing to Two-Stream Instability Observatory

Thank you for your interest in contributing! This document explains the project structure and how to extend it.

---

## Project Structure

```
Two_Stream_Lab/
├── index.html              ← Single-file frontend (theory + live observatory)
├── fetcher.py              ← MMS data downloader and analyser (Python)
├── results/
│   ├── mms_live.json       ← Latest 6-hour window (auto-updated)
│   ├── mms_history.json    ← 6-month rolling record (max 720 entries)
│   └── mms_events.json     ← Unusual event log (max 500 entries)
├── .github/workflows/
│   └── analyze.yml         ← GitHub Actions scheduler (every 6 hours)
├── paper.md                ← JOSE submission paper
├── paper.bib               ← References bibliography
├── CITATION.cff            ← Software citation metadata
└── CONTRIBUTING.md         ← This file
```

---

## The Two Main Components

### 1. `index.html` — Interactive Frontend

The entire frontend is a single self-contained HTML file using:
- **Plotly.js 2.27** (CDN) for all interactive plots
- **Vanilla JavaScript** — no build step, no npm, no bundler
- **CSS custom properties** for the CRT/oscilloscope design system

**To add a new theory section:**
1. Add a `<div class="sec" id="your-section">` block following the existing section template
2. Add a nav link `<a class="nav-link" href="#your-section">Label</a>` in the `<nav>` element
3. Add your JavaScript logic inside the `<script>` block at the bottom

**To update SMILEI PIC data:**
The simulation data is embedded directly in the JavaScript as `const simData = {...}` around line 750. Replace or extend the `v0` keys with new simulation results.

**Key JavaScript functions:**
| Function | Purpose |
|---|---|
| `computeGamma(type, v0, wp, k)` | Core: evaluates Eq.18–21 |
| `genCompare()` | Regenerates the 4-curve comparison plot |
| `plotSimData()` | Overlays SMILEI PIC data vs analytical |
| `drawBox()` | Recomputes SMILEI simulation box parameters |
| `initPhaseSpace()` | Sets up the phase-space canvas animation |

### 2. `fetcher.py` — Data Pipeline

Pure Python (stdlib + `numpy` + `requests`). No other dependencies.

**To add a new instrument or analysis:**
1. Write a `fetch_xxx(t0, t1)` function following the pattern of `fetch_fgm()` — it should return a dict or `None`
2. Call it inside `run()` and add the result to the `live` dict
3. Add corresponding display logic in `index.html` (section 08+)

**To change the fetch cadence:**
Edit `.github/workflows/analyze.yml` — the `cron: '0 */6 * * *'` line sets the 6-hour interval.

---

## Running Locally

**Frontend only (no data):**
```bash
# Just open index.html in any browser — no server needed
open index.html
```

**With live data:**
```bash
pip install numpy requests
python fetcher.py
# Then open index.html — it reads results/*.json from the same directory
```

**GitHub Actions (automated):**
Fork the repo, enable GitHub Pages (Settings → Pages → Deploy from branch → main), and trigger the workflow manually once from the Actions tab. Subsequent runs happen automatically every 6 hours.

---

## Reporting Issues

Please open a GitHub Issue with:
- Which section of the tool is affected (section number and name)
- Browser and OS
- Console errors (F12 → Console tab)
- For `fetcher.py` issues: the full terminal output

---

## Physics Notes for Contributors

- All equations use **normalised units**: $c = 1$, $\omega_{pe} = 1$
- Beam velocity $v_0$ is in units of $c$ (so $v_0 = 0.6$ means $0.6c$)
- The instability is active when the dispersion bracket $[2v_0^2k^2 + \omega_p^2 - \sqrt{8v_0^2k^2\omega_p^2 + \omega_p^4}] < 0$
- MMS data has ~90 day latency; `fetcher.py` automatically targets 90 days ago
- For questions about the physics, open an Issue tagged `physics-question`
