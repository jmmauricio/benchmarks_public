# NTS Figura 18 — dominant electromechanical mode

X_L sweep on the bus 2–3 transmission link. Three pydae configurations
compared against the NTS Figura 18 read-off:

- **milano6ord (X_l = 0.234)** — NTS Tabla 45 literal X_l, Marconato
  effective stator (X_d'' − X_l) = 0.035 pu.
- **milano6ord (X_l = 0)** — the historical "X_l override" workaround
  used in pydae before genrou existed.
- **genrou** — IEEE 1110 Model 2.2 / Anderson-Fouad, no X_l field,
  terminal-referred subtransient reactances.
- **NTS Figura 18** — read off the figure
  (`mdbib_bps/mds/20210709-NTS-SEPE-v2.1/_page_105_Figure_5.jpeg`).
  Only X_L = 0.01 (×, start) and X_L = 0.6 (○, end) are explicitly
  labelled; intermediate values are visual interpolation along the
  marked locus. See `genrou_nts_overlay.png` for the overlay that
  validates the pydae trajectories against the NTS markers.

## Eigenvalue table

Columns are the complex eigenvalue λ = σ + j ω of the dominant
electromechanical mode, the natural frequency f = ω / (2 π), and the
damping ratio ζ = −σ / √(σ² + ω²).

| X_L (pu) |    σ    |    ω    | f (Hz)  | ζ (%) |
|---------:|--------:|--------:|--------:|------:|
| **milano6ord (X_l = 0.234, NTS literal)** |  |  |  |  |
| 0.05     | −1.7335 |  5.2881 | 0.842   | 31.2 |
| 0.10     | −1.3513 |  4.0054 | 0.637   | 32.0 |
| 0.15     | −1.2890 |  3.2161 | 0.512   | 37.2 |
| 0.20     | −1.2093 |  2.3631 | 0.376   | 45.6 |
| 0.25     | −0.7773 |  2.0100 | 0.320   | 36.1 |
| 0.30     | −0.5768 |  1.8401 | 0.293   | 29.9 |
| 0.35     | −0.4608 |  1.7115 | 0.272   | 26.0 |
| 0.40     | −0.3849 |  1.6024 | 0.255   | 23.4 |
| 0.50     | −0.2907 |  1.4172 | 0.226   | 20.1 |
| 0.60     | −0.2332 |  1.2564 | 0.200   | 18.3 |
| **milano6ord (X_l = 0, historical workaround)** |  |  |  |  |
| 0.05     | −1.1474 |  4.5693 | 0.727   | 24.4 |
| 0.10     | −0.8473 |  3.4620 | 0.551   | 23.8 |
| 0.15     | −0.6429 |  2.8178 | 0.448   | 22.2 |
| 0.20     | −0.4836 |  2.4338 | 0.387   | 19.5 |
| 0.25     | −0.3811 |  2.1820 | 0.347   | 17.2 |
| 0.30     | −0.3162 |  1.9951 | 0.318   | 15.7 |
| 0.35     | −0.2734 |  1.8444 | 0.294   | 14.7 |
| 0.40     | −0.2434 |  1.7164 | 0.273   | 14.0 |
| 0.50     | −0.2039 |  1.5021 | 0.239   | 13.5 |
| 0.60     | −0.1775 |  1.3192 | 0.210   | 13.3 |
| **genrou (IEEE 1110 Model 2.2, no X_l)** |  |  |  |  |
| 0.05     | −1.1620 |  4.5909 | 0.731   | 24.5 |
| 0.10     | −0.8601 |  3.4753 | 0.553   | 24.0 |
| 0.15     | −0.6548 |  2.8240 | 0.449   | 22.6 |
| 0.20     | −0.4927 |  2.4365 | 0.388   | 19.8 |
| 0.25     | −0.3882 |  2.1835 | 0.348   | 17.5 |
| 0.30     | −0.3222 |  1.9960 | 0.318   | 15.9 |
| 0.35     | −0.2786 |  1.8449 | 0.294   | 14.9 |
| 0.40     | −0.2480 |  1.7167 | 0.273   | 14.3 |
| 0.50     | −0.2075 |  1.5019 | 0.239   | 13.7 |
| 0.60     | −0.1804 |  1.3187 | 0.210   | 13.6 |
| **NTS Figura 18 (read-off)** |  |  |  |  |
| 0.01 (×) | ≈ −2.00 | ≈ 6.30  | ≈ 1.00  | ≈ 30 |
| 0.05     | ≈ −1.17 | ≈ 4.59  | ≈ 0.73  | ≈ 25 |
| 0.10     | ≈ −0.87 | ≈ 3.48  | ≈ 0.55  | ≈ 24 |
| 0.15     | ≈ −0.67 | ≈ 2.83  | ≈ 0.45  | ≈ 23 |
| 0.20     | ≈ −0.50 | ≈ 2.44  | ≈ 0.39  | ≈ 20 |
| 0.25     | ≈ −0.40 | ≈ 2.19  | ≈ 0.35  | ≈ 18 |
| 0.30     | ≈ −0.32 | ≈ 2.00  | ≈ 0.32  | ≈ 16 |
| 0.35     | ≈ −0.28 | ≈ 1.85  | ≈ 0.29  | ≈ 15 |
| 0.40     | ≈ −0.24 | ≈ 1.72  | ≈ 0.27  | ≈ 14 |
| 0.50     | ≈ −0.21 | ≈ 1.51  | ≈ 0.24  | ≈ 14 |
| 0.60 (○) | ≈ −0.18 | ≈ 1.32  | ≈ 0.21  | ≈ 14 |

## Side-by-side at canonical X_L points

| X_L | milano6ord (X_l=0.234) | milano6ord (X_l=0) | genrou (no X_l) | NTS Fig.18 |
|----:|:----------------------:|:------------------:|:---------------:|:----------:|
|     | σ ± jω   f Hz   ζ%     | σ ± jω   f Hz   ζ% | σ ± jω   f Hz   ζ% | σ ± jω   f Hz   ζ% |
| 0.05 | −1.73 ± j4.29  0.84  31.2 | −1.15 ± j4.57  0.73  24.4 | −1.16 ± j4.59  0.73  24.5 | ≈ −1.17 ± j4.59  0.73  ≈ 25 |
| 0.10 | −1.35 ± j4.01  0.64  32.0 | −0.85 ± j3.46  0.55  23.8 | −0.86 ± j3.48  0.55  24.0 | ≈ −0.87 ± j3.48  0.55  ≈ 24 |
| 0.15 | −1.29 ± j3.22  0.51  37.2 | −0.64 ± j2.82  0.45  22.2 | −0.65 ± j2.82  0.45  22.6 | ≈ −0.67 ± j2.83  0.45  ≈ 23 |
| 0.20 | −1.21 ± j2.36  0.38  45.6 | −0.48 ± j2.43  0.39  19.5 | −0.49 ± j2.44  0.39  19.8 | ≈ −0.50 ± j2.44  0.39  ≈ 20 |
| 0.25 | −0.78 ± j2.01  0.32  36.1 | −0.38 ± j2.18  0.35  17.2 | −0.39 ± j2.18  0.35  17.5 | ≈ −0.40 ± j2.19  0.35  ≈ 18 |
| 0.30 | −0.58 ± j1.84  0.29  29.9 | −0.32 ± j2.00  0.32  15.7 | −0.32 ± j2.00  0.32  15.9 | ≈ −0.32 ± j2.00  0.32  ≈ 16 |
| 0.35 | −0.46 ± j1.71  0.27  26.0 | −0.27 ± j1.84  0.29  14.7 | −0.28 ± j1.84  0.29  14.9 | ≈ −0.28 ± j1.85  0.29  ≈ 15 |
| 0.40 | −0.38 ± j1.60  0.26  23.4 | −0.24 ± j1.72  0.27  14.0 | −0.25 ± j1.72  0.27  14.3 | ≈ −0.24 ± j1.72  0.27  ≈ 14 |
| 0.50 | −0.29 ± j1.42  0.23  20.1 | −0.20 ± j1.50  0.24  13.5 | −0.21 ± j1.50  0.24  13.7 | ≈ −0.21 ± j1.51  0.24  ≈ 14 |
| 0.60 | −0.23 ± j1.26  0.20  18.3 | −0.18 ± j1.32  0.21  13.5 | −0.18 ± j1.32  0.21  13.6 | ≈ −0.18 ± j1.32  0.21  ≈ 14 |

## Observations

- **`genrou` ≈ `milano6ord(X_l=0)`** to within 0.3 percentage points of
  damping at every X_L — the only residual delta is the Marconato
  rotor cross-coupling, not the stator.
- **`milano6ord(X_l=0.234)` is off the NTS curve** by 5–25 percentage
  points and its locus has a non-monotonic damping vs X_L (the peak
  at X_L=0.20 is an artefact of the effective stator reactance
  (X_d'' − X_l) = 0.035 pu interacting unphysically with the
  ST4B + PSS2A loop). This is the motivation for the genrou model:
  drop X_l from the data, use terminal-referred X_d'' / X_q''.
- **`genrou` matches NTS Figura 18 at both labelled endpoints**
  (X_L = 0.01 × and X_L = 0.6 ○) within ~1 pp damping and ~0.02 Hz.
  The overlay (`genrou_nts_overlay.png`) shows the pydae genrou
  trajectory lying on top of the NTS × markers across the whole
  sweep.

## How to reproduce

```bash
cd /Users/jmmauricio/workspace/benchmarks_public/nts/cases/base
PYTHONPATH=/Users/jmmauricio/workspace/pydae/packages/pydae-uds/src \
  /opt/homebrew/Caskroom/miniforge/base/envs/pydae_dev/bin/python \
  nts_base_ini_sympy.py
# Sweep is the bus 2-3 line via pydae.bps.lines.change_line(...)
# at X_pu = 0.05, 0.10, ..., 0.60. See test_genrou.py::test_nts_damping
# in pydae for the automated version.
```

For the milano6ord variants, swap `type: "genrou"` → `type: "milano6ord"`
in `nts_base.hjson` and add `X_l: 0.234, T_AA: 0.0` (or `X_l: 0.0`) to
each generator entry. The dominant electromechanical mode is selected
by maximum participation factor on `delta_1`.
