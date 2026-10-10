# SA-VSS-LMS: Sparse-Aware Adaptive Channel Estimation

A MATLAB simulation framework for estimating sparse multipath channels with adaptive filters. It implements and fairly compares zero-attracting, proportionate, and robust LMS-family algorithms, and introduces **CZA-LMS (Cosine Zero-Attracting LMS)**, a zero-attraction method with a raised-cosine attraction function.

The long-term goal is **SA-VSS-LMS**: an algorithm that uses an online estimate of channel sparsity to adapt both the step size μ and the attraction strength ρ.

> **Status:** CZA-LMS is a proposed method. Its novelty has not been confirmed by a complete literature review, and all results below come from simulation only.

---

## Problem setup

A known training sequence `d[n]` passes through an unknown sparse FIR channel `h` (most taps are zero) and is received with noise:

```text
x[n] = sum_k h_k[n] d[n-k] + v[n]
```

The adaptive filter is configured for **channel identification**:

```text
adaptive-filter input = d[n]      (known training sequence)
desired response      = x[n]      (received signal)
estimated weights     = channel estimate
```

So the weights estimate the channel directly, and the main metric compares them with the true channel:

$$\text{NMSD}(n) = 10\log_{10}\frac{\lVert \mathbf{w}(n)-\mathbf{h}(n)\rVert^2}{\lVert \mathbf{h}(n)\rVert^2}\ \text{dB}$$

(This is different from equalisation, where the filter would estimate an inverse of the channel and NMSD against `h` would not be meaningful.)

---

## Algorithms

Every algorithm uses the same update structure:

$$\mathbf{w}(n+1) = \mathbf{w}(n) + \mu\, f\big(e(n)\big)\,\mathbf{x}(n) - \rho\, c\big(\mathbf{w}(n)\big)$$

- **f(e)**, the error (score) function, sets how the filter reacts to the current error. Changing it gives robustness to outliers.
- **c(w)**, the attraction function (applied per tap), pulls taps toward zero. Changing it gives sparsity awareness.

| Family | Algorithm | File | f(e) or gain | c(w) |
|---|---|---|---|---|
| Baseline | LMS | `lms_c.m` | e | — |
| Baseline | NLMS | `nlms.m` | e / ‖x‖² | — |
| Baseline | VSS-LMS | `vss_lms.m` | e, with μ(n+1) = αμ(n) + γe² | — |
| Zero-attracting | ZA-LMS | `za_lms.m` | e | sign(w) |
| Zero-attracting | RZA-LMS | `rza_lms.m` | e | sign(w) / (1 + ε\|w\|) |
| Zero-attracting | l0-LMS | `l0_lms.m` | e | sign(w)(1 − β\|w\|) for \|w\| ≤ 1/β |
| Zero-attracting | PZA-LMS | `pza_lms.m` | e | sign(w)[1 − (α−1)\|w\|] / (1+\|w\|)^(α+1) for \|w\| ≤ 1/(α−1) |
| Zero-attracting | **CZA-LMS (proposed)** | `cza_lms.m` | e | sign(w)·½[1 + cos(π\|w\|/a)] for \|w\| ≤ a |
| Proportionate | PNLMS | `pnlms.m` | per-tap gain ∝ \|w_l\| | — |
| Proportionate | IPNLMS | `ipnlms.m` | blend of uniform and proportionate gains | — |
| Robust | MVC (Versoria) | `mvc.m` | pτ\|e\|^(p−1)sign(e) / (1 + τ\|e\|^p)² | — |
| Robust | Andrew's sine | `andrew_sine.m` | sin(e/c) for \|e\| ≤ πc, else 0 | — |
| Robust | HCAF | `hcaf.m` | sinh(λe) | — |
| Robust | LHCAF | `lhcaf.m` | tanh(λe) | — |

### CZA-LMS

$$c(w) = \begin{cases}\operatorname{sign}(w)\,\tfrac{1}{2}\left[1+\cos\!\left(\dfrac{\pi |w|}{a}\right)\right], & |w|\le a\\[4pt] 0, & |w|>a\end{cases}$$

Design properties, which follow from the formula:

- **`a` sets only the width.** Writing u = |w|/a gives ½(1 + cos πu) for every a, so curves for different a are identical up to scaling. PZA's α changes both the cutoff and the shape.
- **Flat top.** Taps just off zero get nearly the full pull.
- **Smooth landing.** Both the value and the slope are zero at |w| = a, so taps near the cutoff are not hit by an abrupt change in pull.
- **Behaves like a smoothed hard threshold.** It pulls harder than l0-LMS's linear taper in the inner half of the zone and softer in the outer half. The implied penalty saturates at a/2, like a truncated l0 count.
- **Cost:** O(M) per iteration, with one comparison per tap and one cosine per tap inside the zone. This is cheaper than PZA's high power, but more expensive than l0-LMS's single multiply.

---

## Repository structure

```text
main.m                   Single run of all configured algorithms; defines cfg and alg
runEqualizer.m           Dispatches an algorithm by name
evaluateResult.m         NMSD, MSD and MSE trajectories and final values
generateSparseChannel.m  Random sparse channel generator
cza_lms.m                Proposed CZA-LMS algorithm

equalizers/              All other adaptive filters (see table above)
environments/            Input signal, channel models, noise, received-data generation

rho_sweep.m              Tunes rho for each zero-attracting algorithm
monte_carlo_study.m      Main comparison over hundreds of random sparse channels
small_tap_study.m        How each attraction shape treats real taps near the cutoff
ablation_study.m         Effect of CZA's cutoff a for user-defined channels
```

### Environment factors

The channel and receiver are built from interchangeable parts set in `cfg`:

```matlab
cfg.signal.generator    = @whiteGaussianSignal;
cfg.channel.initializer = @sparseChannel;
cfg.channel.factors     = {};                 % e.g. {@slowDriftChannel} or {@abruptChangeChannel}
cfg.noise.apply         = @awgnNoise;
```

A static sparse channel is a sparse initializer with no factors, and a drifting one adds `@slowDriftChannel`. Combinations need no new environment files.

---

## Quick start

1. Open MATLAB in the repository root. The scripts add subfolders to the path automatically.
2. Run **`main.m`** for a single run of all configured algorithms, with NMSD, MSD and MSE curves.
3. Run **`monte_carlo_study.m`** for the full comparison. It runs `rho_sweep.m` first if `rho_sweep_results.mat` is missing.

Each study script calls `main.m` once to load `cfg` and `alg`. The block at the end of the parameters section in `main.m` makes it return early in that case:

```matlab
if isappdata(0, 'SAVSS_INIT_ONLY') && getappdata(0, 'SAVSS_INIT_ONLY')
    return;
end
```

All study scripts are seeded, so results are reproducible.

---

## Study scripts

### `rho_sweep.m`

Runs ZA, RZA, l0, PZA and CZA at ρ ∈ {10⁻⁵, 3×10⁻⁵, 10⁻⁴, 3×10⁻⁴, 10⁻³} on 200 random channels and picks the ρ with the lowest mean final NMSD for each. Results are saved to `rho_sweep_results.mat`. It uses different channels from the Monte Carlo study, so ρ is not tuned on the test set. A warning is printed if a best ρ lands on the edge of the grid.

### `monte_carlo_study.m`

Compares all algorithms on 400 random sparse channels, using each sparse algorithm's tuned ρ.

- **Channels:** K active taps (K uniform from 2 to 25% of M), random positions, magnitudes uniform in [0.2, 1], random signs.
- **Pairing:** every algorithm sees the same channel, input and noise in each run.
- **Outputs:** mean final NMSD with 95% confidence interval, spread across channels, convergence times, number of diverged runs, a per-K breakdown, and averaged NMSD curves.

Set `forceRetune = true` at the top to rerun the ρ sweep. The script warns if the saved ρ values were tuned under different settings.

### `small_tap_study.m`

Places one real tap of magnitude s inside or near the shared cutoff (0.1) and sweeps s. For LMS, l0-LMS, PZA and CZA it reports:

- final NMSD, and the bias on the small tap compared with the theoretical prediction;
- a **thresholding view**: how much of a far-below-cutoff tap is removed, how much of a near-cutoff tap is kept, and leftover energy on true-zero taps.

`st.rhoMode = "common"` uses the same ρ for every method and isolates the effect of shape. `"tuned"` uses each method's best ρ.

### `ablation_study.m`

Sweeps CZA's cutoff a for channels defined directly in the script:

```matlab
chans(end+1) = struct('taps', [2 4 9], 'gains', [0.2 0.8 0.7], 'label', "");
```

It reports the safe window for a (within 1 dB of the best NMSD) for each channel, and plots NMSD against a/|h|min to test whether the window scales with the smallest true tap.

---

## Fair-comparison setup

- **Matched effective step size.** Near convergence, every robust update reduces to (μ × shape factor) × e. Each μ is set so that this effective step equals 0.01, e.g. μ = 0.01/λ for LHCAF and μ = 0.01·c for Andrew's sine. NLMS, PNLMS and IPNLMS use μ = 0.01·M·σx² for the same average step per tap.
- **Matched cutoffs.** CZA's a = 0.1, l0-LMS's 1/β = 0.1 and PZA's 1/(α−1) = 0.1, so their attraction zones are identical.
- **Tuned ρ** for each zero-attracting method, chosen on held-out channels.
- **Averaging.** NMSD is averaged in linear scale over channels, then converted to dB.

### Metrics

| Metric | Meaning |
|---|---|
| `FinalNMSD_dB` | Mean NMSD over the last 10% of iterations, averaged over channels |
| `CI95` | 95% confidence interval of that mean. Overlapping intervals mean the difference is not established |
| `Std_dB` | Spread of final NMSD across channels |
| `T_target` | First iteration at which the averaged NMSD curve reaches −30 dB |
| `T_conv` | First iteration within 3 dB of the algorithm's own final level |
| `Diverged` | Runs that diverged, excluded from averages |

---

## Current findings

Static channel, white Gaussian input, AWGN at 30 dB SNR, M = 32, N = 20,000 iterations. Rerun the scripts to regenerate the exact numbers.

- **Setup check.** With matched effective steps, LMS, NLMS, MVC, Andrew's sine, HCAF and LHCAF all reach about −37 dB, as expected in Gaussian noise. Their robustness only matters under impulsive noise.
- **Zero attraction.** Over 400 random channels, **CZA-LMS, PZA-LMS and l0-LMS all reach about −45 dB, and their 95% confidence intervals overlap**, about 8 dB better than LMS. RZA reaches about −42 dB and ZA about −39 dB. All sparse gains shrink as the channel becomes less sparse.
- **Proportionate.** PNLMS and IPNLMS converge quickly but settle at higher error floors than NLMS, the known trade-off of proportionate gains.
- **Cutoff matters more than shape.** In `small_tap_study.m`, choosing a cutoff above a real tap costs up to about 10 dB, while the difference between attraction shapes is about 2–4 dB.
- **Thresholding behaviour.** With tuned ρ, CZA removes more of far-below-cutoff taps **and** keeps more of near-cutoff taps than l0-LMS. PZA keeps near-cutoff taps best. Removing a real tap is not free: its energy appears as extra noise on the other taps.
- **Practical rule for a.** Choose a below the smallest tap considered significant and above the noise-level jitter of zero taps. Inside that window, final NMSD is insensitive to the exact value.

---

## Roadmap

- [ ] SA-VSS-LMS: drive μ(n) and ρ(n) from an online sparsity estimate (e.g. the Hoyer measure)
- [ ] Time-varying channels (`slowDriftChannel`, `abruptChangeChannel`) in the Monte Carlo study
- [ ] Impulsive noise model (e.g. Bernoulli–Gaussian or α-stable)
- [ ] Combined robust and sparse filters, e.g. CZA attraction with an LHCAF error function
- [ ] Literature review to confirm the novelty of CZA-LMS
- [ ] Runtime and operation-count comparison

---

## References

- Y. Chen, Y. Gu, A. O. Hero, "Sparse LMS for system identification," *ICASSP*, 2009. (ZA-LMS, RZA-LMS)
- Y. Gu, J. Jin, S. Mei, "l0 norm constraint LMS algorithm for sparse system identification," *IEEE Signal Processing Letters*, 2009.
- D. L. Duttweiler, "Proportionate normalized least-mean-squares adaptation in echo cancelers," *IEEE Trans. Speech and Audio Processing*, 2000.
- J. Benesty, S. L. Gay, "An improved PNLMS algorithm," *ICASSP*, 2002.
- F. Huang, J. Zhang, S. Zhang, "Maximum Versoria criterion-based robust adaptive filtering algorithm," *IEEE Trans. Circuits and Systems II*, 2017.
- S. Wang, W. Wang, K. Xiong, H. H. C. Iu, C. K. Tse, "Logarithmic hyperbolic cosine adaptive filter and its performance analysis," *IEEE Access*, 2019.
- PZA-LMS: Maheshwari and George, IIT Gandhinagar, 2016 (add full citation).