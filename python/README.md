# KPC-Toolbox (Python)

Python edition of the KPC-Toolbox for fitting Markovian Arrival Processes (MAPs)
and phase-type (PH) distributions to empirical traces via Kronecker Product
Composition (KPC).

This is a thin, standalone package: it does **not** reimplement any algorithm
that already ships in [`line-solver`](https://pypi.org/project/line-solver/).
Every fitting routine and numerical primitive is imported from `line_solver` and
re-exported under a flat namespace that mirrors the MATLAB KPC-Toolbox API
(`kpcfit_*`, `map_*`, `trace_*`, `aph_*`, `mmpp2_*`, `ctmc_*`/`dtmc_*`). On top of
that, this package adds the material the MATLAB toolbox ships that has no home in
`line-solver`: the example MAP library, runnable demos with bundled traces, and
tests.

## Install

```bash
pip install kpc-toolbox        # pulls in line-solver, numpy, scipy
```

## Quick start

```python
import numpy as np
import kpctoolbox as k

# --- Phase-type fitting: match the first moments of a distribution ---
D0, D1 = k.map_hyperexp(1.0, 8.0, 0.6)          # a hyper-exponential, SCV=8
E = [float(k.map_moment(D0, D1, i)) for i in (1, 2, 3)]
options = k.kpcfit_ph_options(np.array(E))
PH = k.kpcfit_ph_auto(np.array(E), options)      # list of fitted PH (D0, D1)

# --- MAP fitting: fit autocorrelation of a trace ---
from kpctoolbox import demo
S = demo.load_trace("BCAUG89")                   # interarrival-time samples
trace = k.kpcfit_init(S)
MAP = k.kpcfit_auto(trace, OnlyAC=True)          # fitted MAP {D0, D1}
```

## What's here

- `kpctoolbox` — facade re-exporting the KPC engine and MAP/trace/PH primitives
  from `line-solver`.
- `kpctoolbox.examples` — the `maplib_*` example MAP library (Erlang,
  hyper-exponential, MMPP, circulant, saw-tooth, and the pre-fitted Bellcore
  Aug89 MAP).
- `kpctoolbox.demo` — runnable demos (`demo_kpcfit_bcaug89`, `demo_kpcfit_ph`)
  and trace loaders (`load_trace`), with the BCAUG89, DEC-PKT-1-UDP, and
  LIVEMAPS traces bundled as compressed `.npz`.

## Entry points

| Task           | Functions                                        |
|----------------|--------------------------------------------------|
| MAP fitting    | `kpcfit_init` → `kpcfit_auto`                     |
| PH fitting     | `kpcfit_ph_options` → `kpcfit_ph_auto`           |

Note: line-solver's PH-fitting options use snake_case keyword arguments
(e.g. `kpcfit_ph_options(E, max_num_states=8)`), unlike the MATLAB PascalCase.

## References

1. G. Casale, E. Z. Zhang, E. Smirni. *KPC-Toolbox: Best Recipes for Automatic
   Trace Fitting Using Markovian Arrival Processes.* Performance Evaluation,
   67(9):873-896, 2010.
2. G. Casale, E. Z. Zhang, E. Smirni. *Trace Data Characterization and Fitting
   for Markov Modeling.* Performance Evaluation, 67(2):61-79, 2010.

Released under the BSD-3 license (see `LICENSE.TXT`).
