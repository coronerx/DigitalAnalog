# Lab 7 & 8: Forced Oscillations — Analysis Code

This repository contains the MATLAB analysis scripts used to produce the results and
figures in the *Lab 7 & 8: Forced Oscillations* writeup. This README summarizes the
experiment, the data pipeline, and how the two scripts relate to each other.

## Files

| File | Description |
|---|---|
| `disc_analysis.m` | Analysis for the **large damping disc** configuration, 5 motor voltages |
| `plate_analysis.m` | Analysis for the **small paper plate** configuration, 11 motor voltages, plus an undamped "mass only" trial |
| `partA*.txt`, `partB1-*.txt` | Lab 7 data: free (`A`) and driven (`B`) oscillation trials (not included here — see original data folder) |
| `partf*.txt`, `partg*.txt` | Lab 8 data: free (`f1`–`f3`), driven (`f11`–`f1_11`), and undamped driven (`g1`–`g11`) trials |

Each raw `.txt` data file is logged from an ultrasonic position sensor via Logger Pro,
with columns `[t, x, v, a, ?, omega, alpha]` — time, position, velocity, acceleration,
and (for driven trials) the motor's angular velocity/acceleration.

## Experimental setup

A mass hangs from a spring, driven at one end by a motor-and-arm mechanism that can
pull and release the spring at a controllable rate (see Figure 1 in the writeup). An
ultrasonic position sensor records the mass's vertical position over time. Two damping
elements were tested in separate labs:

- **Lab 7:** a large damping disc (150 g)
- **Lab 8:** a small paper plate (10.40 g)

For each configuration, the scripts analyze:
1. **Free oscillation** (motor off) to measure the natural frequency ω₀ and compare it
   to the theoretical value ω₀ = √(k/m).
2. **Driven oscillation** (motor on, several fixed voltages) to fit the damping
   resistance, amplitude, phase lag, and quality factor as functions of driving
   frequency ω.

## Physics / model

The system is modeled as a driven, damped harmonic oscillator:

```
m·a = F0·sin(ωt) − k·x − R·v        (Eq. 1)
```

where `m` is the oscillating mass, `k` the spring constant, `F0` the driving force
amplitude, `R` the velocity-proportional damping resistance, and `ω` the driving
angular frequency. Solving for `R` at each data sample gives:

```
R = |m·a − F0·sin(ωt) + k·x| / v
```

The damping rate is `γ = R / (2m)`. Given `γ` and the measured `ω₀`, the steady-state
response is predicted by the standard driven-oscillator solution:

```
A(ω)   = (F0/m) / sqrt[(ω0² − ω²)² + 4γ²ω²]      — amplitude
φ(ω)   = atan[2γω / (ω² − ω0²)]                   — phase lag
v_max  = ω·A(ω)                                    — velocity extrema
```

The quality factor is estimated from the half-power bandwidth (the range of ω over
which the response stays above ~70% of its peak):

```
Q = ω0 / (2·Δω)                                     (Eq. 2)
```

Uncertainty in `R` (and hence `γ` and `Q`) is propagated from the uncertainties in
`F0`, `ω`, `a`, `k`, `x`, and `v` via standard error propagation (Eq. 3 in the
writeup), since `R` is a function of all of these measured quantities.

## Analysis pipeline (both scripts)

1. **Load data** for each trial into per-column vectors/matrices (one column per
   voltage for the driven trials).
2. **Estimate ω from free oscillation:** find velocity zero-crossings (turning points
   of x), keep every other crossing to isolate full periods, average the time between
   crossings 2 apart to get the period `T`, then `ω = 2π/T`.
3. **Fit damping:** apply Eq. 1 pointwise across all samples to solve for `R`, then
   average (with outlier/Inf rejection where `v ≈ 0`) to get `γ` per voltage.
4. **Compute predicted A(ω), φ(ω), v_max(ω)** and compare to the measured data
   (Figures 2–4 for Lab 7 disc; Figures 8–11 for Lab 8 plate).
5. **Estimate Q** from the half-power bandwidth of `v_max` vs. `ω`.
6. **Propagate uncertainty** through `R`, `γ`, and `Q`.
7. **Normalize:** rescale `ω` by `ω0` and plot the dimensionless response curves
   `f(ω/ω0)`, `φ(ω/ω0)`, and `g(ω/ω0)` so that the disc and plate configurations can
   be compared directly on the same axes (Figures 15–17).

## Key differences between the two scripts

- **Voltage count:** Lab 7 handles 5 voltage trials with individually-named
  variables (`tB11`, `tB12`, …); Lab 8 handles 11 trials loaded in a `for` loop into
  matrices, since hand-naming 11 sets of variables would be unwieldy.
- **Bandwidth search:** Lab 7 finds the single closest sample to the 70%-of-peak
  velocity level on each side of resonance. Lab 8 uses a more robust version that
  collects *all* candidate crossings within a tolerance and picks the closest
  rising/falling pair — this matters more with 11 trials, which have noisier,
  less resonance-like curves (see writeup: none of the Lab 8 driven trials show a
  clear turning point in `v_max`, meaning the system is over/critically damped at all
  11 tested frequencies).
- **Error propagation loop:** Lab 7 computes the six error terms of Eq. 3 with
  vectorized `mean()` calls. Lab 8 does the equivalent computation with explicit
  nested loops and Inf-rejection, since some `v_C` samples are exactly zero.
- **Part G (Lab 8 only):** an additional undamped driven-oscillation trial
  (`partg*.txt`), used to show the amplitude growing over time near resonance
  (Figure 14), since with negligible damping the motor continuously pumps energy
  into the system.

## Known quirks / caveats (kept intentional, not "bugs" to fix)

These are called out in code comments rather than corrected, since correcting them
would change the numbers reported in the writeup:

- In Lab 7's Part-A section, the A2 and A3 blocks concatenate trial vectors with
  commas (row-wise) while A1 uses semicolons (column-wise) — likely a copy/paste
  inconsistency, but left as-is.
- `m` (total oscillating mass) is reassigned partway through each script to a
  hand-picked value (`m=0.51`) that differs slightly from `m_d+m_a` — this is used
  only for the final normalized-curve section.
- Both scripts print some intermediate values (`Rerr`, `gammaErr`, `gamma`) to the
  console via `disp`/`fprintf` for manual inspection while iterating on the analysis.

## Summary of reported results

| Quantity | Lab 7 (disc) | Lab 8 (plate) |
|---|---|---|
| Damping mass | 150 g | 10.40 g |
| ω₀ (measured, rad/s) | 5.946 | 5.915 |
| F₀ (N) | 0.917 | 1.480 |
| γ (s⁻¹) | 11.679 | 8.606 |

The disc produces significantly stronger damping (larger γ) than the plate, as
expected given its larger surface area / mass. Across both configurations, all
computed γ were at or above ω₀ (i.e. the system was heavily/over-damped), which is
why the amplitude and velocity response curves show a monotonic decline rather than a
clear resonance peak, and why measured driving frequencies never reached the
theoretical ω₀.

## Requirements

- MATLAB (tested informally; no toolboxes beyond base MATLAB are required)
- Raw `.txt` data files in the same working directory as the scripts
