"""
Two-Stream Instability Observatory — MMS Data Fetcher
=======================================================
Vivek Shrivastav, Sikkim University (2025)

Downloads real MMS spacecraft data via NASA HAPI and:
  1. Fetches the latest 6-hour window of L2 survey data
  2. Computes all plasma parameters relevant to TSI detection
  3. Validates measured values against theoretical predictions (Eq.18-21)
  4. Detects unusual events and electron holes
  5. Maintains a 6-month rolling history

NASA HAPI endpoint: https://cdaweb.gsfc.nasa.gov/hapi
MMS L2 data latency: ~90 days (we fetch 90 days ago to guarantee availability)

Instruments:
  FGM  — Magnetic field B (GSE), 8 sps survey, dataset MMS1_FGM_SRVY_L2@0
  EDP  — Electric field E (GSE), 32 sps, dataset MMS1_EDP_FAST_L2_DCE
  FPI-DES — Electron moments (~4.5 s), dataset MMS1_FPI_FAST_L2_DES-MOMS
  FPI-DIS — Ion moments (~4.5 s),      dataset MMS1_FPI_FAST_L2_DIS-MOMS
"""

import json, os, traceback, math
import numpy as np
import requests
from datetime import datetime, timezone, timedelta

LIVE_FILE    = "results/mms_live.json"
HISTORY_FILE = "results/mms_history.json"
EVENTS_FILE  = "results/mms_events.json"
HAPI         = "https://cdaweb.gsfc.nasa.gov/hapi"
SC           = "1"

# ── Time window ────────────────────────────────────────────────────────────────
def get_trange():
    """6-hour window, 90 days ago (guaranteed to be available)."""
    now   = datetime.now(timezone.utc)
    stop  = now - timedelta(days=90)
    start = stop - timedelta(hours=6)
    fmt   = "%Y-%m-%dT%H:%M:%SZ"
    return start.strftime(fmt), stop.strftime(fmt), start, stop

# ── HAPI helpers ───────────────────────────────────────────────────────────────
def _hapi(dataset, parameter, t0, t1):
    url = f"{HAPI}/data"
    p   = {"id":dataset,"parameters":parameter,"time.min":t0,"time.max":t1,"format":"csv"}
    try:
        r = requests.get(url, params=p, timeout=90)
        if r.status_code != 200:
            print(f"    HTTP {r.status_code} — {dataset}/{parameter}")
            return None
        lines = [l for l in r.text.strip().split('\n') if not l.startswith('#')]
        rows  = []
        for line in lines:
            parts = line.split(',')
            try:
                vals = [float(v) for v in parts[1:] if v.strip()]
                if vals: rows.append(vals)
            except: pass
        if rows: print(f"    {dataset}/{parameter}: {len(rows)} rows")
        return rows if rows else None
    except Exception as e:
        print(f"    {dataset}/{parameter} error: {e}"); return None

def col(rows, i=0):
    """Extract column, remove fill values (|val|>1e29)."""
    if not rows: return np.array([])
    a = np.array([r[i] for r in rows if i < len(r)], dtype=float)
    a[np.abs(a) > 1e29] = np.nan
    return a[np.isfinite(a)]

def sub(a, n=600):
    if len(a) == 0: return []
    s = max(1, len(a)//n)
    return [round(float(v),5) for v in np.asarray(a)[::s]]

# ── FGM — Magnetic field ───────────────────────────────────────────────────────
def fetch_fgm(t0, t1):
    print("\n  [FGM]")
    # Query catalog first to find working dataset name
    rows = None
    try:
        cat = requests.get(f"{HAPI}/catalog", timeout=30).json()
        fgm_ids = [d["id"] for d in cat.get("catalog",[]) if "MMS1" in d["id"] and "FGM" in d["id"]]
    except:
        fgm_ids = []

    candidates = [(i, "mms1_fgm_b_gse_srvy_l2_clean") for i in fgm_ids]
    candidates += [
        ("MMS1_FGM_SRVY_L2@0", "mms1_fgm_b_gse_srvy_l2_clean"),
        ("MMS1_FGM_SRVY_L2",   "mms1_fgm_b_gse_srvy_l2_clean"),
        ("MMS1_FGM_SRVY_L2@0", "mms1_fgm_b_gse_srvy_l2"),
    ]
    for ds, param in candidates:
        rows = _hapi(ds, param, t0, t1)
        if rows:
            print(f"    Using {ds}/{param}")
            break

    if not rows: return None

    nc = len(rows[0])
    if nc >= 4:
        Bx=col(rows,0); By=col(rows,1); Bz=col(rows,2); Bm=col(rows,3)
    elif nc >= 3:
        Bx=col(rows,0); By=col(rows,1); Bz=col(rows,2)
        n=min(len(Bx),len(By),len(Bz))
        Bm=np.sqrt(Bx[:n]**2+By[:n]**2+Bz[:n]**2)
    elif nc >= 1:
        Bm=col(rows,0); Bx=By=Bz=np.zeros_like(Bm)
    else: return None

    ok = (Bm > 0) & (Bm < 1e5) & np.isfinite(Bm)
    if ok.sum() < 20: return None
    Bm_ok = Bm[ok]

    # Turbulence PSD
    N=len(Bm_ok); dt=0.125  # FGM survey: 8 sps
    win=np.hanning(N)
    Bf=np.fft.rfft((Bm_ok-Bm_ok.mean())*win)
    freqs=np.fft.rfftfreq(N,d=dt)
    PSD=np.abs(Bf)**2/(np.sum(win**2)*dt)
    ok_f=(freqs>0.01)&(freqs<2.0)&(PSD>0)&np.isfinite(PSD)
    alpha=None
    if ok_f.sum()>=8:
        c,_=np.polyfit(np.log(freqs[ok_f]),np.log(PSD[ok_f]),1)
        alpha=round(float(c),3)

    print(f"    OK: B={Bm_ok.mean():.2f}±{Bm_ok.std():.2f} nT  α={alpha}")
    return {
        "n_points":       int(ok.sum()),
        "B_mean_nT":      round(float(Bm_ok.mean()),3),
        "B_std_nT":       round(float(Bm_ok.std()),3),
        "B_max_nT":       round(float(Bm_ok.max()),3),
        "B_min_nT":       round(float(Bm_ok.min()),3),
        "spectral_index": alpha,
        "turbulence_regime": (
            "Kolmogorov-like (α ≈ −5/3)"  if alpha and abs(alpha+5/3)<0.2
            else "Iroshnikov-Kraichnan (α ≈ −3/2)" if alpha and abs(alpha+1.5)<0.2
            else f"Steeper kinetic range (α={alpha})" if alpha and alpha < -2.0
            else f"Other (α={alpha})" if alpha else "Not determined"
        ),
        "Bmag":  sub(Bm_ok),
        "Bx":    sub(Bx[ok] if len(Bx)==len(Bm) else np.array([])),
        "By":    sub(By[ok] if len(By)==len(Bm) else np.array([])),
        "Bz":    sub(Bz[ok] if len(Bz)==len(Bm) else np.array([])),
        "psd_freq": [round(float(v),6) for v in freqs[1:300]],
        "psd_vals": [round(float(v),8) for v in PSD[1:300]],
    }

# ── EDP — Electric field ───────────────────────────────────────────────────────
def fetch_edp(t0, t1):
    print("\n  [EDP]")
    rows = _hapi("MMS1_EDP_FAST_L2_DCE","mms1_edp_dce_gse_fast_l2",t0,t1)
    if not rows: return None
    nc=len(rows[0])
    if nc>=3:
        Ex=col(rows,0); Ey=col(rows,1); Ez=col(rows,2)
    elif nc>=1:
        Ez=col(rows,0); Ex=Ey=np.zeros_like(Ez)
    else: return None
    n=min(len(Ex),len(Ey),len(Ez))
    Emag=np.sqrt(Ex[:n]**2+Ey[:n]**2+Ez[:n]**2); Epar=Ez[:n].copy()
    ok=(Emag>0)&(Emag<1e5)&np.isfinite(Emag)
    if ok.sum()<50: return None
    Em=Emag[ok]; Ep=Epar[ok]; dt=1/32.0

    # ── Electron hole detection ──────────────────────────────────────────────
    Ep_abs=np.abs(Ep); valid=Ep_abs[Ep_abs>0]
    if len(valid)==0: return None
    thr=max(3.0*float(np.median(valid)),0.5)
    n_holes=0; hole_events=[]; fs=1/dt; i=1
    while i<len(Ep)-1:
        if Ep[i]>thr:
            j=i
            while j>max(0,i-int(0.05*fs)) and Ep[j]>0: j-=1
            k2=i
            while k2<min(len(Ep)-1,i+int(0.05*fs)) and Ep[k2]>0: k2+=1
            if k2>j and Ep[k2]<=0:
                dur=(k2-j)*dt*1e3
                if 0.1<dur<500:
                    n_holes+=1
                    if len(hole_events)<20:
                        hole_events.append({"E_max_mVm":round(float(Ep[i]),3),
                                            "duration_ms":round(dur,2)})
                    i=k2+1; continue
        i+=1

    # E-field PSD
    N_f=len(Em); win_f=np.hanning(N_f)
    Efft=np.fft.rfft((Em-Em.mean())*win_f)
    fE=np.fft.rfftfreq(N_f,d=dt)
    PSDE=np.abs(Efft)**2/(np.sum(win_f**2)*dt)

    print(f"    OK: E={Em.mean():.3f} mV/m  holes={n_holes}  thr={thr:.3f}")
    return {
        "n_points":    int(ok.sum()),
        "E_mean_mVm":  round(float(Em.mean()),4),
        "E_std_mVm":   round(float(Em.std()),4),
        "E_max_mVm":   round(float(Em.max()),4),
        "Epar_std_mVm":round(float(np.std(Ep)),4),
        "bipolar_threshold_mVm": round(thr,3),
        "electron_holes": {
            "n_detected": n_holes,
            "status": (
                "LIKELY ELECTRON HOLES" if n_holes>5
                else "POSSIBLE HOLES" if n_holes>0
                else "NONE DETECTED"
            ),
            "note": (
                f"{n_holes} bipolar E∥ events above {thr:.2f} mV/m — "
                +("consistent with nonlinear TSI saturation products (electron holes/phase-space vortices)"
                  if n_holes>5
                  else "marginal — may be noise or weak TSI activity"
                  if n_holes>0
                  else "no clear bipolar signatures in this window")
            ),
            "events": hole_events
        },
        "Emag_series": sub(Em),
        "Epar_series": sub(Ep),
        "langmuir_psd": {
            "freq_Hz":     [round(float(v),4) for v in fE[1:300]],
            "PSD_mVm2_Hz": [round(float(v),8) for v in PSDE[1:300]]
        }
    }

# ── FPI-DES — Electron moments ─────────────────────────────────────────────────
def fetch_fpi_electrons(t0, t1):
    print("\n  [FPI-DES]")
    result={}; DS="MMS1_FPI_FAST_L2_DES-MOMS"

    rn=_hapi(DS,"mms1_des_numberdensity_fast",t0,t1)
    if rn:
        Ne=col(rn,0); ok=(Ne>0)&(Ne<1e4)
        if ok.sum()>5:
            Ne_ok=Ne[ok]
            fpe=8980.0*np.sqrt(Ne_ok.mean())
            result.update({
                "Ne_mean_cc":  round(float(Ne_ok.mean()),4),
                "Ne_std_cc":   round(float(Ne_ok.std()),4),
                "Ne_max_cc":   round(float(Ne_ok.max()),4),
                "Ne_series":   sub(Ne_ok),
                "fpe_Hz":      round(float(fpe),1),
                "fpe_kHz":     round(float(fpe/1e3),3),
                # Debye length λ_D = v_th / ω_pe
                # Will be computed after temperature is known
            })

    rv=_hapi(DS,"mms1_des_bulkv_gse_fast",t0,t1)
    if rv and len(rv[0])>=3:
        Vx=col(rv,0); Vy=col(rv,1); Vz=col(rv,2)
        n=min(len(Vx),len(Vy),len(Vz))
        Vm=np.sqrt(Vx[:n]**2+Vy[:n]**2+Vz[:n]**2)
        ok=(Vm>0)&(Vm<5e4)&np.isfinite(Vm)
        if ok.sum()>5:
            Vm_ok=Vm[ok]
            result.update({
                "Ve_mean_kms": round(float(Vm_ok.mean()),2),
                "Ve_std_kms":  round(float(Vm_ok.std()),2),
                "Ve_max_kms":  round(float(Vm_ok.max()),2),
                "Ve_series":   sub(Vm_ok),
                "beam_flag":   bool(Vm_ok.mean()>300),
                "beam_note":   (
                    f"ELECTRON BEAM — mean Ve={Vm_ok.mean():.1f} km/s > 300 km/s threshold. "
                    "TSI conditions present."
                    if Vm_ok.mean()>300
                    else f"No clear beam — mean Ve={Vm_ok.mean():.1f} km/s < 300 km/s threshold."
                )
            })

    rtp=_hapi(DS,"mms1_des_temppara_fast",t0,t1)
    rte=_hapi(DS,"mms1_des_tempperp_fast",t0,t1)
    if rtp and rte:
        Tp=col(rtp,0); Te=col(rte,0); n=min(len(Tp),len(Te))
        if n>5:
            Tp_ok=Tp[:n]; Te_ok=Te[:n]; ok=(Tp_ok>0)&(Te_ok>0)
            if ok.sum()>5:
                Tf=Tp_ok[ok]; Tef=Te_ok[ok]
                ratio=Tef/Tf; A=ratio-1.0
                c_kms=3e5
                # Thermal speed: v_th = sqrt(2 k_B T_par / m_e) in km/s
                # m_e c^2 = 511 keV = 511000 eV → v_th = sqrt(2T[eV]/511000) × c
                vth=np.sqrt(2*float(Tf.mean())/511e3)*c_kms
                Ve_m=result.get("Ve_mean_kms",0)
                dfdv=Ve_m/max(vth,1.0)

                # Debye length: λ_D = v_th / (sqrt(2) × ω_pe) in metres
                # ω_pe = 2π × f_pe
                fpe_hz=result.get("fpe_Hz",1e4)
                omega_pe=2*np.pi*fpe_hz
                vth_ms=vth*1e3
                debye_m=vth_ms/(np.sqrt(2)*omega_pe) if omega_pe>0 else None

                result.update({
                    "Tpar_mean_eV":   round(float(Tf.mean()),2),
                    "Tperp_mean_eV":  round(float(Tef.mean()),2),
                    "Tratio_mean":    round(float(ratio.mean()),4),
                    "anisotropy_A":   round(float(A.mean()),4),
                    "Tpar_series":    sub(Tf),
                    "Tperp_series":   sub(Tef),
                    "Tratio_series":  sub(ratio),
                    "vth_par_kms":    round(vth,2),
                    "dfdv_proxy":     round(dfdv,4),
                    "debye_length_m": round(debye_m,4) if debye_m else None,
                    "dfdv_note": (
                        f"v_beam/v_th = {dfdv:.2f} >> 1 → BUMP-ON-TAIL present, df/dv > 0 at v_phase. "
                        "Classical two-stream mechanism active."
                        if dfdv > 3
                        else f"v_beam/v_th = {dfdv:.2f} — marginal. Possible truncated distribution (Jebaraj+2025)."
                        if dfdv > 1.5
                        else f"v_beam/v_th = {dfdv:.2f} < 1.5 → no bump-on-tail. "
                        "Beam insufficient for classical TSI."
                    )
                })

                # Weibel susceptibility
                A_mean=float(A.mean())
                if A_mean>0:
                    vp=np.sqrt(2*float(Tef.mean())/511e3)*c_kms
                    g_wp=np.sqrt(A_mean/(1+A_mean))*vp/c_kms
                    result["weibel_gamma_over_wp"]=round(g_wp,6)
                result["weibel_susceptibility"]=(
                    f"HIGH — A={A_mean:.3f} >> 0. Weibel instability strongly driven." if A_mean>0.3
                    else f"MODERATE — A={A_mean:.3f}. Marginal Weibel conditions." if A_mean>0.05
                    else f"LOW — A={A_mean:.3f} ≈ 0. Plasma near isotropic." if A_mean>-0.05
                    else f"FIREHOSE — A={A_mean:.3f} < 0. Firehose instability possible."
                )

    print(f"    OK: Ne={result.get('Ne_mean_cc','?')} cc  "
          f"Ve={result.get('Ve_mean_kms','?')} km/s  "
          f"Tpar={result.get('Tpar_mean_eV','?')} eV  "
          f"Tratio={result.get('Tratio_mean','?')}")
    return result if result else None

# ── FPI-DIS — Ion moments ──────────────────────────────────────────────────────
def fetch_fpi_ions(t0, t1):
    print("\n  [FPI-DIS]")
    result={}; DS="MMS1_FPI_FAST_L2_DIS-MOMS"
    rn=_hapi(DS,"mms1_dis_numberdensity_fast",t0,t1)
    if rn:
        Ni=col(rn,0); ok=(Ni>0)&(Ni<1e4)
        if ok.sum()>5: result["Ni_mean_cc"]=round(float(Ni[ok].mean()),4)
    rtp=_hapi(DS,"mms1_dis_tempperp_fast",t0,t1)
    rta=_hapi(DS,"mms1_dis_temppara_fast",t0,t1)
    if rtp and rta:
        Tp=col(rtp,0); Ta=col(rta,0); n=min(len(Tp),len(Ta))
        if n>5:
            ok=(Tp[:n]>0)&(Ta[:n]>0)
            if ok.sum()>5:
                result.update({
                    "Ti_perp_eV": round(float(Tp[:n][ok].mean()),2),
                    "Ti_par_eV":  round(float(Ta[:n][ok].mean()),2),
                    "Ti_ratio":   round(float((Tp[:n][ok]/Ta[:n][ok]).mean()),4)
                })
    if result: print(f"    OK: Ni={result.get('Ni_mean_cc','?')} cc")
    return result if result else None

# ── Theory Validation ──────────────────────────────────────────────────────────
def validate_theory(fgm, fpi_e, fpi_i):
    """
    Compare measured MMS plasma conditions against theoretical TSI predictions.

    Uses equations from Shrivastav et al. (2025):
      Eq.(18) Non-relativistic single-species:
        γ = √{ ½[2v₀²k² + ω²ₚ − √(8v₀²k²ω²ₚ + ω⁴ₚ)] }
      Eq.(19) Relativistic single-species:
        γ = √{ ½[2v₀²k² + ω²ₚ/γ₀³ − √(8v₀²k²ω²ₚ/γ₀³ + ω⁴ₚ/γ₀⁶)] }

    v₀ = V_e / c  (electron beam velocity in units of c)
    ωₚ = 1        (normalised plasma frequency = 1 by definition)
    k  = k_res = ωₚ / v₀ = 1/v₀  (resonance condition k·v₀ = ωₚ)
    """
    if not fpi_e:
        return {"available": False, "note": "No FPI-DES data for theory validation"}

    Ve_kms = fpi_e.get("Ve_mean_kms")
    vth_kms= fpi_e.get("vth_par_kms")
    Ne_cc  = fpi_e.get("Ne_mean_cc")
    fpe_Hz = fpi_e.get("fpe_Hz")

    if not all([Ve_kms, Ne_cc, fpe_Hz]):
        return {"available": False, "note": "Insufficient parameters for theory validation"}

    c_kms  = 3e5
    v0     = Ve_kms / c_kms           # beam velocity normalised to c
    wp     = 1.0                       # normalised plasma frequency
    k_res  = wp / v0 if v0 > 0 else 0 # resonant wavenumber

    # Lorentz factor
    if v0 >= 1.0:
        return {"available": False, "note": f"v₀={v0:.4f} ≥ c — relativistic overflow"}
    g0 = 1.0 / np.sqrt(1 - v0**2)

    results = {}

    # Sweep k around resonance: k_res ± 50%
    k_min = max(0.001, k_res * 0.1)
    k_max = k_res * 3.0 if k_res > 0 else 5.0
    N_k   = 300
    k_arr = np.linspace(k_min, k_max, N_k)

    for eq_name, multi, relativistic in [
        ("eq18_nr_single", False, False),
        ("eq19_rel_single",False, True),
        ("eq20_nr_multi",  True,  False),
        ("eq21_rel_multi", True,  True),
    ]:
        gr_arr = []
        for k in k_arr:
            v2=v0**2; k2=k**2; wp2=wp**2; wp4=wp2**2
            if relativistic:
                g3=g0**3; g6=g0**6
                if multi:
                    t1=2*v2*k2+2*wp2/g3; inner=16*v2*k2*wp2/g3+4*wp4/g6
                else:
                    t1=2*v2*k2+wp2/g3;   inner=8*v2*k2*wp2/g3+wp4/g6
            else:
                if multi:
                    t1=2*v2*k2+2*wp2; inner=16*v2*k2*wp2+4*wp4
                else:
                    t1=2*v2*k2+wp2;   inner=8*v2*k2*wp2+wp4
            bracket = t1 - np.sqrt(inner)
            gr = np.sqrt(0.5*abs(bracket)) if bracket < 0 else 0.0
            gr_arr.append(round(float(gr),6))

        max_gr  = max(gr_arr) if gr_arr else 0
        k_at_max= float(k_arr[np.argmax(gr_arr)]) if gr_arr else 0

        # Real growth rate in physical units: γ_phys = γ_norm × ω_pe
        gamma_phys_rad_s = max_gr * 2 * np.pi * fpe_Hz

        results[eq_name] = {
            "max_growth_rate_normalised": round(max_gr, 6),
            "k_at_max":                  round(k_at_max, 4),
            "max_growth_rate_rad_s":     round(gamma_phys_rad_s, 2),
            # e-folding time: τ = 1/γ_phys in seconds
            "efolding_time_ms":          round(1e3/gamma_phys_rad_s, 4) if gamma_phys_rad_s > 0 else None,
            "unstable_range_k":          [round(k_min,4), round(k_max,4)],
            "growth_curve_k":            [round(float(k),4) for k in k_arr[::3]],
            "growth_curve_gr":           gr_arr[::3],
        }

    # Plasma conditions summary
    # β = plasma beta: ratio thermal to magnetic pressure
    B_nT = fgm.get("B_mean_nT") if fgm else None
    beta = None
    if B_nT and Ne_cc and fpi_e.get("Tpar_mean_eV"):
        beta = round(0.403 * Ne_cc * fpi_e["Tpar_mean_eV"] / B_nT**2, 4)

    # Alfvén speed
    vA = None
    if B_nT and fpi_i and fpi_i.get("Ni_mean_cc"):
        vA = round(21.8 * B_nT / np.sqrt(fpi_i["Ni_mean_cc"]), 2)

    # Is TSI active? Check if measured conditions satisfy instability criterion
    tsi_active = bool(v0 > 0.001 and k_res > 0 and results.get("eq18_nr_single",{}).get("max_growth_rate_normalised",0) > 0)

    return {
        "available":        True,
        "input_parameters": {
            "v0_over_c":   round(v0, 6),
            "lorentz_g0":  round(float(g0), 6),
            "k_resonance": round(k_res, 4),
            "omega_pe_rad_s": round(2*np.pi*fpe_Hz, 2),
            "fpe_kHz":     fpi_e.get("fpe_kHz"),
            "Ne_cc":       Ne_cc,
            "Ve_kms":      Ve_kms,
            "B_nT":        B_nT,
            "plasma_beta": beta,
            "alfven_speed_kms": vA,
        },
        "tsi_active_predicted": tsi_active,
        "theory_curves":   results,
        "validation_note": (
            f"Theory predicts TSI active for v₀/c={v0:.4f}: "
            f"max γ = {results['eq18_nr_single']['max_growth_rate_normalised']:.5f} ωₚ "
            f"at k={results['eq18_nr_single']['k_at_max']:.3f}. "
            f"E-field doubles every {results['eq18_nr_single']['efolding_time_ms']:.2f} ms."
            if tsi_active and results.get("eq18_nr_single",{}).get("max_growth_rate_normalised",0)>0
            else "Beam velocity too low for classical TSI at observed conditions."
        )
    }

# ── Langmuir proxy ─────────────────────────────────────────────────────────────
def langmuir_proxy(edp, fpi_e):
    if not edp or not fpi_e:
        return {"active": False, "note": "No data"}
    fpe=fpi_e.get("fpe_Hz")
    psd_d=edp.get("langmuir_psd",{})
    freqs=np.array(psd_d.get("freq_Hz",[])); vals=np.array(psd_d.get("PSD_mVm2_Hz",[]))
    if fpe is None or len(freqs)<5:
        return {"active":False,"note":"f_pe unavailable","fpe_Hz":fpe}
    band=(freqs>0.8*fpe)&(freqs<1.2*fpe)&(vals>0)
    bkg =(freqs>0.5*fpe)&(freqs<0.7*fpe)&(vals>0)
    if band.sum()<2 or bkg.sum()<2:
        return {"active":False,"note":f"f_pe={fpe/1e3:.2f} kHz outside EDP frequency range (Nyquist ~16 Hz for slow survey)","fpe_Hz":fpe}
    rat=float(np.mean(vals[band]))/max(float(np.mean(vals[bkg])),1e-30)
    act=bool(rat>3.0)
    return {
        "active": act,
        "fpe_Hz": round(fpe,1), "fpe_kHz": round(fpe/1e3,3),
        "E_at_fpe":     round(float(np.mean(vals[band])),8),
        "E_background": round(float(np.mean(vals[bkg])),8),
        "enhancement_ratio": round(rat,3),
        "note": (
            f"Langmuir waves ACTIVE — E-field enhanced {rat:.1f}× near f_pe={fpe/1e3:.2f} kHz. "
            "Consistent with beam-driven Langmuir wave generation."
            if act else
            f"No significant enhancement near f_pe={fpe/1e3:.2f} kHz (R={rat:.2f} < 3)."
        )
    }

# ── TSI→Weibel chain ───────────────────────────────────────────────────────────
def detect_chain(fpi_e):
    if not fpi_e: return {"detected":False,"note":"No data"}
    Ve_s=fpi_e.get("Ve_series",[]); Tr_s=fpi_e.get("Tratio_series",[])
    A=fpi_e.get("anisotropy_A",0)
    if len(Ve_s)<10 or len(Tr_s)<10:
        return {"detected":False,"note":"Insufficient time series"}
    n=min(len(Ve_s),len(Tr_s))
    Ve=np.array(Ve_s[:n]); Tr=np.array(Tr_s[:n])
    Ve_n=(Ve-Ve.mean())/(Ve.std()+1e-10); Tr_n=(Tr-Tr.mean())/(Tr.std()+1e-10)
    ml=min(n//4,20)
    xc=[float(np.correlate(Ve_n[:n-l],Tr_n[l:])[0]/(n-l)) for l in range(ml)]
    pl=int(np.argmax(xc)); pc=xc[pl]
    det=(fpi_e.get("beam_flag",False) and A>0.05 and pc>0.2 and pl>0)
    return {
        "detected": det,
        "peak_correlation": round(pc,4),
        "peak_lag_samples": pl,
        "beam_flag": fpi_e.get("beam_flag",False),
        "anisotropy_A": A,
        "note": (
            f"TSI→Weibel chain DETECTED — V_e spike leads T⊥/T∥ rise by {pl} samples "
            f"(cross-correlation r={pc:.3f} > 0.2). "
            "Consistent with PIC simulation predictions (Innocenti et al. 2017)."
            if det else
            "Chain not detected — " + (
                "no electron beam present" if not fpi_e.get("beam_flag")
                else f"anisotropy too low (A={A:.3f} < 0.05)" if A<=0.05
                else f"cross-correlation weak (r={pc:.3f} < 0.2)"
            )
        ),
        "xcorr": [round(v,4) for v in xc]
    }

# ── Unusual event detection ────────────────────────────────────────────────────
def detect_unusual_events(fgm, edp, fpi_e, fpi_i, theory, trange_start):
    """
    Detect and classify unusual plasma events in this 6-hour window.
    Events are flagged and stored for the events log.
    """
    events = []

    # 1. Large magnetic field (magnetospheric compression/shock)
    if fgm and fgm.get("B_max_nT") and fgm.get("B_mean_nT"):
        B_max = fgm["B_max_nT"]
        B_mean= fgm["B_mean_nT"]
        B_std = fgm.get("B_std_nT",0)
        # Flag if max > mean + 3σ, or if absolute value very high (> 500 nT)
        if B_max > B_mean + 3*B_std or B_max > 500:
            events.append({
                "type": "HIGH_B_FIELD",
                "severity": "HIGH" if B_max > 500 else "MODERATE",
                "time": trange_start,
                "values": {"B_max_nT": B_max, "B_mean_nT": B_mean},
                "explanation": (
                    f"Unusually large magnetic field: B_max = {B_max:.1f} nT "
                    f"({(B_max-B_mean)/max(B_std,0.1):.1f}σ above mean). "
                    "Possible causes: magnetospheric compression, flux tube boundary, "
                    "current sheet crossing, or reconnection inflow region."
                )
            })

    # 2. Many electron holes
    n_holes = (edp or {}).get("electron_holes",{}).get("n_detected",0)
    if n_holes >= 10:
        events.append({
            "type": "ELECTRON_HOLE_BURST",
            "severity": "HIGH" if n_holes >= 20 else "MODERATE",
            "time": trange_start,
            "values": {"n_holes": n_holes},
            "explanation": (
                f"{n_holes} electron holes (bipolar E∥ events) detected in 6 hours. "
                "This is unusually high activity, indicating an active or recently saturated "
                "two-stream instability. These phase-space vortices are the direct nonlinear "
                "saturation products of Langmuir wave trapping (Bernstein et al. 1957, "
                "Graham et al. 2016 GRL)."
            )
        })

    # 3. High temperature anisotropy
    A = (fpi_e or {}).get("anisotropy_A", 0)
    if A > 0.3:
        events.append({
            "type": "HIGH_TEMPERATURE_ANISOTROPY",
            "severity": "HIGH" if A > 0.5 else "MODERATE",
            "time": trange_start,
            "values": {"A": round(A,4), "Tratio": (fpi_e or {}).get("Tratio_mean")},
            "explanation": (
                f"High electron temperature anisotropy: A = T⊥/T∥ − 1 = {A:.3f}. "
                "This drives the Weibel instability — electrons preferentially heated "
                "in the perpendicular direction spontaneously form current filaments "
                "that generate magnetic fields. This can follow TSI saturation when "
                "Langmuir waves scatter beam electrons into perpendicular velocities."
            )
        })

    # 4. High electron bulk speed (beam)
    Ve = (fpi_e or {}).get("Ve_mean_kms", 0)
    if Ve > 300:
        events.append({
            "type": "ELECTRON_BEAM_DETECTED",
            "severity": "HIGH" if Ve > 600 else "MODERATE",
            "time": trange_start,
            "values": {"Ve_mean_kms": Ve, "Ve_max_kms": (fpi_e or {}).get("Ve_max_kms")},
            "explanation": (
                f"Fast electron beam: mean V_e = {Ve:.1f} km/s > 300 km/s threshold. "
                f"Ratio v_beam/v_th = {(fpi_e or {}).get('dfdv_proxy','?')}. "
                "When v_beam >> v_thermal, the electron velocity distribution function "
                "develops a positive slope df/dv > 0 at v_phase = ω_pe/k, "
                "which is the classical two-stream instability condition. "
                "Theory (Eq.18) predicts maximum growth rate γ_max at k_res = ω_pe/v_beam."
            )
        })

    # 5. TSI→Weibel chain detected
    if (theory or {}).get("tsi_active_predicted") and A > 0.05:
        events.append({
            "type": "TSI_WEIBEL_CHAIN",
            "severity": "MODERATE",
            "time": trange_start,
            "values": {
                "v0_over_c": (theory.get("input_parameters",{}) or {}).get("v0_over_c"),
                "A": round(A,4),
                "predicted_gamma": (theory.get("theory_curves",{}).get("eq18_nr_single",{}) or {}).get("max_growth_rate_normalised")
            },
            "explanation": (
                "Both TSI and Weibel conditions are simultaneously active. "
                "PIC simulations (Innocenti et al. 2017) show: beam drives Langmuir waves "
                "→ Langmuir waves scatter electrons into ⊥ direction → builds T⊥ > T∥ "
                "→ drives Weibel electromagnetic waves. This is one of the first observational "
                "windows where this chain may be directly tested with MMS data."
            )
        })

    # 6. Extreme turbulence
    alpha = (fgm or {}).get("spectral_index")
    if alpha and alpha < -2.5:
        events.append({
            "type": "STEEP_TURBULENCE_SPECTRUM",
            "severity": "MODERATE",
            "time": trange_start,
            "values": {"spectral_index": alpha},
            "explanation": (
                f"Unusually steep magnetic power spectrum: α = {alpha:.3f} "
                f"(much steeper than Kolmogorov α = -5/3 = -1.667). "
                "This indicates enhanced dissipation at sub-ion scales, possibly due to "
                "Landau damping of kinetic Alfvén waves or reconnection-driven turbulence. "
                "Active Langmuir/TSI waves can inject energy into the perpendicular "
                "direction, feeding into the turbulent cascade."
            )
        })

    return events

# ── Main ───────────────────────────────────────────────────────────────────────
def run():
    os.makedirs("results", exist_ok=True)
    now = datetime.now(timezone.utc)
    t0, t1, dt_start, dt_stop = get_trange()
    errors = []

    print("="*65)
    print("TWO-STREAM INSTABILITY OBSERVATORY — MMS DATA FETCHER")
    print(f"Data window : {t0}  →  {t1}")
    print(f"Spacecraft  : MMS-{SC}")
    print(f"API         : {HAPI}")
    print("="*65)

    fgm=edp=fpi_e=fpi_i=None

    for label, fn in [
        ("FGM",     lambda: fetch_fgm(t0,t1)),
        ("EDP",     lambda: fetch_edp(t0,t1)),
        ("FPI-DES", lambda: fetch_fpi_electrons(t0,t1)),
        ("FPI-DIS", lambda: fetch_fpi_ions(t0,t1)),
    ]:
        try:
            r=fn()
            if label=="FGM":     fgm=r
            elif label=="EDP":   edp=r
            elif label=="FPI-DES": fpi_e=r
            elif label=="FPI-DIS": fpi_i=r
        except Exception as e:
            errors.append(f"{label}: {str(e)}"); print(f"    ERROR: {e}")

    # Derived analyses
    theory   = validate_theory(fgm, fpi_e, fpi_i)
    chain    = detect_chain(fpi_e)
    langmuir = langmuir_proxy(edp, fpi_e)
    events   = detect_unusual_events(fgm, edp, fpi_e, fpi_i, theory, t0)

    # Derived quantities
    derived={}
    if fgm and fpi_i and fpi_i.get("Ni_mean_cc"):
        derived["alfven_speed_kms"]=round(float(21.8*fgm["B_mean_nT"]/np.sqrt(fpi_i["Ni_mean_cc"])),2)
    if fgm and fpi_e and fpi_e.get("Ne_mean_cc") and fpi_e.get("Tpar_mean_eV") and fgm["B_mean_nT"]>0:
        derived["plasma_beta"]=round(0.403*fpi_e["Ne_mean_cc"]*fpi_e["Tpar_mean_eV"]/fgm["B_mean_nT"]**2,4)

    # Overall TSI status
    n_holes = (edp or {}).get("electron_holes",{}).get("n_detected",0)
    dfdv    = (fpi_e or {}).get("dfdv_proxy",0)
    A       = (fpi_e or {}).get("anisotropy_A",0)
    beam    = (fpi_e or {}).get("beam_flag",False)

    tsi_score = 0
    if beam:                                                    tsi_score += 2
    if n_holes > 5:                                             tsi_score += 3
    if n_holes > 0:                                             tsi_score += 1
    if langmuir.get("active"):                                  tsi_score += 2
    if chain.get("detected"):                                   tsi_score += 3
    if A > 0.05:                                                tsi_score += 1
    if dfdv > 3:                                                tsi_score += 2

    tsi_status = (
        "ACTIVE"    if tsi_score >= 7
        else "PROBABLE" if tsi_score >= 4
        else "POSSIBLE" if tsi_score >= 2
        else "QUIET"
    )

    live = {
        "fetched_utc":      now.isoformat(),
        "data_window_utc":  [t0, t1],
        "data_window_note": "MMS L2 data ~90 days behind real-time",
        "spacecraft":       f"MMS-{SC}",
        "tsi_status":       tsi_status,
        "tsi_score":        tsi_score,
        "unusual_events":   events,
        "magnetic_field":   fgm,
        "electric_field":   edp,
        "electrons":        fpi_e,
        "ions":             fpi_i,
        "theory_validation":theory,
        "tsi_weibel_chain": chain,
        "langmuir_proxy":   langmuir,
        "derived":          derived,
        "errors":           errors
    }
    with open(LIVE_FILE,"w") as f: json.dump(live,f,indent=2)
    print(f"\nSaved → {LIVE_FILE}")

    # ── History ────────────────────────────────────────────────────────────────
    try:
        with open(HISTORY_FILE) as f: hist=json.load(f)
        recs=hist.get("records",[])
    except FileNotFoundError:
        recs=[]

    recs.append({
        "utc":           now.isoformat(),
        "data_utc":      t0,
        "B_mean_nT":     fgm.get("B_mean_nT")       if fgm   else None,
        "B_alpha":       fgm.get("spectral_index")   if fgm   else None,
        "Ne_cc":         fpi_e.get("Ne_mean_cc")     if fpi_e else None,
        "fpe_kHz":       fpi_e.get("fpe_kHz")        if fpi_e else None,
        "Ve_kms":        fpi_e.get("Ve_mean_kms")    if fpi_e else None,
        "Tpar_eV":       fpi_e.get("Tpar_mean_eV")   if fpi_e else None,
        "Tperp_eV":      fpi_e.get("Tperp_mean_eV")  if fpi_e else None,
        "Tratio":        fpi_e.get("Tratio_mean")    if fpi_e else None,
        "A":             fpi_e.get("anisotropy_A")   if fpi_e else None,
        "dfdv":          fpi_e.get("dfdv_proxy")     if fpi_e else None,
        "n_holes":       n_holes if edp else None,
        "beam":          beam,
        "langmuir":      langmuir.get("active"),
        "chain":         chain.get("detected"),
        "tsi_status":    tsi_status,
        "tsi_score":     tsi_score,
        "n_events":      len(events),
    })
    recs=recs[-720:]  # 6 months at 6-hour cadence = 720 records
    with open(HISTORY_FILE,"w") as f:
        json.dump({"last_updated_utc":now.isoformat(),"n_records":len(recs),"records":recs},f,indent=2)
    print(f"History → {len(recs)} records (~{len(recs)*6/24:.0f} days)")

    # ── Events log ─────────────────────────────────────────────────────────────
    try:
        with open(EVENTS_FILE) as f: ev_log=json.load(f)
        all_events=ev_log.get("events",[])
    except FileNotFoundError:
        all_events=[]

    all_events.extend(events)
    all_events=all_events[-500:]  # keep last 500 events
    with open(EVENTS_FILE,"w") as f:
        json.dump({"last_updated_utc":now.isoformat(),"n_events":len(all_events),"events":all_events},f,indent=2)
    print(f"Events  → {len(all_events)} total events logged")

    print(f"\n{'='*65}")
    print(f"TSI STATUS: {tsi_status} (score={tsi_score}/13)")
    print(f"  FGM:      {'OK  B='+str(fgm.get('B_mean_nT'))+'nT  α='+str(fgm.get('spectral_index')) if fgm else 'FAILED'}")
    print(f"  EDP:      {'OK  holes='+str(n_holes) if edp else 'FAILED'}")
    print(f"  FPI-DES:  {'OK  Ne='+str(fpi_e.get('Ne_mean_cc'))+'cc  Ve='+str(fpi_e.get('Ve_mean_kms'))+'km/s' if fpi_e else 'FAILED'}")
    print(f"  FPI-DIS:  {'OK  Ni='+str(fpi_i.get('Ni_mean_cc'))+'cc' if fpi_i else 'FAILED'}")
    print(f"  Beam:     {beam}")
    print(f"  Chain:    {chain.get('detected')}")
    print(f"  Langmuir: {langmuir.get('active')}")
    print(f"  Unusual events this window: {len(events)}")
    if errors: print(f"  Errors: {errors}")

if __name__ == "__main__":
    run()
