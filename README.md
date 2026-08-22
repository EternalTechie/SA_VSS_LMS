# SA-VSS-LMS Simulation Project

This project is a MATLAB simulation framework for studying adaptive filtering over sparse multipath communication channels. The long-term goal is to develop and evaluate a Sparse-Aware Variable Step-Size Least Mean Squares (SA-VSS-LMS) algorithm against baselines such as Normalized Least Mean Squares (NLMS).

The current focus is the simulation environment, not the SA-VSS-LMS algorithm itself.

## Current Setup

The simulator currently models a discrete-time finite impulse response (FIR) multipath channel:

```text
x[n] = sum_k h_k[n] d[n-k] + v[n]
```

where `d[n]` is the known input/training signal, `h[n]` is the true channel impulse response, `v[n]` is additive noise, and `x[n]` is the received signal.

The default input signal is a real-valued white Gaussian discrete-time sequence. Binary Phase Shift Keying (BPSK) is intentionally not used in the current setup because the first objective is adaptive-filter behavior and channel recovery, not digital modulation performance.

## Why Channel Identification

The current NLMS baseline is configured for channel identification:

```text
adaptive-filter input = d[n]
desired response      = x[n]
estimated weights     = channel estimate
```

This was chosen because the main research question is sparse channel recovery. In this configuration, the adaptive weights estimate the true channel directly, so Normalized Mean Square Deviation (NMSD) between the estimated channel and the true channel is meaningful.

This is different from equalization, where the adaptive filter would receive `x[n]` and try to recover `d[n]`. In that case, the weights estimate an inverse filter rather than the channel itself, so directly comparing the weights to `h[n]` would not be a valid channel-recovery metric.

## Modular Architecture

The code is organized so that input generation, channel/receiver simulation, adaptive filters, and evaluation metrics can be changed independently:

```text
input generator -> channel/receiver factors -> data -> adaptive filter -> result -> evaluation
```

The main files are:

- `main.m`: top-level orchestration only.
- `generateInputSignal.m`: generates the known discrete-time input/training sequence.
- `generateReceivedData.m`: produces the received signal from the supplied input and configured channel/receiver factors.
- `generateChannelTrajectory.m`: builds the true channel trajectory from an initializer and optional channel factors.
- `runEqualizer.m`: selects the adaptive filter.
- `sparseChannel.m`: creates the current sparse initial finite impulse response channel.
- `slowDriftChannel.m`: optional channel factor that adds slow random-walk drift to active taps.
- `applyTimeVaryingFIR.m`: applies the true finite impulse response channel trajectory to the input signal.
- `awgnNoise.m`: adds Additive White Gaussian Noise (AWGN).
- `equalizers/`: adaptive filtering algorithms.
- `evaluateResult.m`: common evaluation metrics.
- `metrics/`: reusable metric helpers.

The input signal `d[n]` is generated before the channel/receiver stage. This keeps the excitation signal independent from the simulated environment. The channel/receiver stage receives `d[n]`, generates the true channel trajectory `h[n]`, applies the channel, adds receiver noise, and returns the received signal `x[n]`.

The returned data structure contains the input signal, received signal, true channel trajectory, and any channel-change indices. The true channel is stored because this is a simulation and objective channel-estimation metrics require ground truth. Adaptive algorithms do not receive the true channel.

Each adaptive filter returns a common result structure containing output, error, and coefficient trajectory. This allows future algorithms such as variable-step-size NLMS, zero-attracting NLMS, and SA-VSS-LMS to use the same evaluation path.

## Evaluation

The primary current metric is NMSD:

```text
NMSD[n] = ||w[n] - h[n]||^2 / ||h[n]||^2
```

where `w[n]` is the estimated channel and `h[n]` is the true channel. The evaluator also reports error power as a diagnostic quantity, but error power is not the primary channel-recovery metric.

## Factor Composition

The channel/receiver simulation is moving toward composable factors rather than one function per scenario. The current configuration uses function handles:

```matlab
cfg.signal.generator = @whiteGaussianSignal;

cfg.channel.initializer = @sparseChannel;
cfg.channel.factors = {};              % Example: {@slowDriftChannel}

cfg.noise.apply = @awgnNoise;
```

This means a static sparse channel is not a separate environment function. It is a sparse initial channel with no time-varying channel factors. A sparse channel with slow drift can use the same sparse initializer plus `@slowDriftChannel` in `cfg.channel.factors`.

This avoids creating a new environment file for every combination of sparse structure, slow drift, abrupt changes, anomalies, fractional delay, and noise.

## Design Direction

Channel structure, channel evolution, fractional delay, noise, anomalies, and receiver impairments should continue to be selectable through configuration so that combinations can be tested without duplicating environment code.

This incremental structure is intentional: first verify the simulator and NLMS baseline, then study step-size behavior and metrics, and only later add sparsity-aware and variable-step-size mechanisms.
