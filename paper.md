---
title: 'Two-Stream Instability Observatory: An Interactive Theory Laboratory and Live NASA MMS Data Monitor'
tags:
  - plasma physics
  - two-stream instability
  - Weibel instability
  - PIC simulation
  - SMILEI
  - NASA MMS
  - dispersion relation
  - magnetosphere
  - electron holes
  - JavaScript
  - Python
authors:
  - name: Vivek Shrivastav
    affiliation: 1
affiliations:
  - name: Department of Physics, Sikkim University, Sikkim, India
    index: 1
date: 2025
bibliography: paper.bib
---

# Summary

The two-stream instability (TSI) is a fundamental collective phenomenon in plasma physics, occurring when two counter-propagating electron beams transfer energy to electrostatic waves, leading to exponential growth of Langmuir wave amplitudes. Understanding TSI is critical across astrophysics, magnetospheric physics, laser-plasma interactions, and inertial confinement fusion research [@Buneman1958; @Weibel1959].

The **Two-Stream Instability Observatory** is an open-source, browser-based research platform consisting of two tightly integrated components: a **Theory Laboratory** providing interactive exploration of exact analytical dispersion relations validated against SMILEI particle-in-cell (PIC) simulations [@Derouillat2018], and a **Live MMS Observatory** that downloads real NASA Magnetospheric Multiscale (MMS-1) spacecraft data every six hours and tests whether theoretical TSI conditions are active in the Earth's magnetosphere [@Burch2016].

# Statement of Need

Teaching and researching plasma instabilities presents a persistent challenge: the gap between analytical theory and observable data is wide, the mathematics is non-trivial, and computational tools (PIC codes, data reduction pipelines) require significant setup. No existing tool bridges all three simultaneously — interactive dispersion theory, PIC simulation comparison, and live spacecraft data validation — in a single zero-install environment.

Existing resources each address only part of the problem. Standalone PIC codes such as SMILEI [@Derouillat2018] require a local HPC environment and Python post-processing expertise. Browser-based particle simulations (e.g., particleincell.com) show raw particle dynamics but provide no analytical growth rate predictions or spacecraft data. NASA's CDAWeb interface provides MMS data but requires custom analysis pipelines. This tool closes all three gaps simultaneously for students, educators, and researchers alike.

# Theory Laboratory

The laboratory implements four closed-form dispersion relations derived from the two-fluid equations for counter-propagating symmetric electron beams (Shrivastav et al., 2025, in preparation):

**Equation 18 — Non-relativistic single-species (ions frozen):**
$$\gamma = \sqrt{\frac{1}{2}\left[2v_0^2k^2 + \omega_p^2 - \sqrt{8v_0^2k^2\omega_p^2 + \omega_p^4}\right]}$$

**Equation 19 — Relativistic single-species:**
$$\gamma = \sqrt{\frac{1}{2}\left[2v_0^2k^2 + \frac{\omega_p^2}{\gamma_0^3} - \sqrt{\frac{8v_0^2k^2\omega_p^2}{\gamma_0^3} + \frac{\omega_p^4}{\gamma_0^6}}\right]}$$

**Equations 20–21** extend the above to mobile-ion (multi-species) cases, replacing $\omega_p^2$ with $2\omega_p^2$ in the dispersion bracket. Here $\gamma_0 = 1/\sqrt{1-v_0^2}$ is the Lorentz factor and normalised units ($c = \omega_{pe} = 1$) are used throughout.

All four equations are implemented as real-time JavaScript computations with interactive sliders for beam velocity $v_0$, plasma frequency $\omega_p$, and wavenumber $k$. Growth rate curves $\gamma(k)$ are rendered via Plotly [@Plotly]. Equations 18–21 are validated against embedded SMILEI PIC simulation data for $v_0 = 0.1c$ to $0.95c$, showing agreement within 5% across all regimes. The laboratory includes a Simulation Box Designer that computes SMILEI-compatible parameters (box length, cell resolution, CFL condition, saturation time estimate) from user inputs.

# MMS Live Observatory

A Python pipeline (`fetcher.py`, 814 lines) runs automatically every six hours via GitHub Actions, downloading MMS-1 Level 2 data from the NASA CDAS HAPI API [@Annex2024]. Four instruments are queried:

| Instrument | Dataset | Quantity | Cadence |
|---|---|---|---|
| FGM | `MMS1_FGM_SRVY_L2` | Magnetic field **B** (GSE) | 8 sps |
| EDP | `MMS1_EDP_FAST_L2_DCE` | Electric field **E** (GSE) | 32 sps |
| FPI-DES | `MMS1_FPI_FAST_L2_DES-MOMS` | $n_e$, $V_e$, $T_\parallel$, $T_\perp$ | ~4.5 s |
| FPI-DIS | `MMS1_FPI_FAST_L2_DIS-MOMS` | $n_i$, $T_i$ | ~4.5 s |

From these measurements, the pipeline computes: (1) theoretical TSI growth rates using Equations 18–21 with the measured $v_0 = V_e/c$ and $\omega_{pe} = \sqrt{n_e e^2/m_e\epsilon_0}$; (2) electron hole detection via bipolar $E_\parallel$ pulse identification [@Graham2016; @Steinvall2021]; (3) bump-on-tail proxy $v_\text{beam}/v_\text{thermal}$ [@Jebaraj2025]; (4) TSI→Weibel causal chain test via cross-correlation of $V_e(t)$ and $T_\perp/T_\parallel(t)$ [@Innocenti2017]; (5) Langmuir wave proxy from E-field PSD near $f_{pe}$; and (6) magnetic turbulence spectral index $\alpha$ from FGM power spectra [@Kolmogorov1941].

Results are stored as rolling JSON files (up to 720 records, approximately six months of history) committed automatically to the repository, and rendered live by the dashboard at [vivek-shrivastav.github.io/Two_Stream_Lab](https://vivek-shrivastav.github.io/Two_Stream_Lab/).

# Educational Design

The tool serves both research and pedagogy. Section 01 includes a Beginner/Expert mode toggle, presenting the same physics at two levels of mathematical depth. FAQ entries with expandable answers cover nine conceptual questions spanning instability onset, the role of relativity, multi-species dynamics, SMILEI parameter mapping, and the distinction between TSI and Weibel instabilities. A Parameter Playground allows users to save, compare, and overlay multiple simulation runs with CSV export. Phase-space animation visualises particle trapping and cat's-eye vortex formation in real time. The six-month historical record and unusual event log expose students to real observational data patterns across different magnetospheric regions.

# Acknowledgements

The author thanks the SMILEI development team at Ecole Polytechnique for the open-source PIC code, and NASA/GSFC for public access to MMS Level 2 data via the CDAS HAPI interface. Supervisory guidance from the Department of Physics, Sikkim University is gratefully acknowledged.

# References
