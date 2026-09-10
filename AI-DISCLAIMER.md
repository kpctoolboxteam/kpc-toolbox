# AI Disclaimer

The Python package under `python/` (PyPI `kpc-toolbox`, `import kpctoolbox`) is a
port of the MATLAB toolbox produced with Claude Opus 5 (Anthropic), used through
Claude Code. This note states that plainly, so that users can weigh it when
deciding how far to rely on the code.

## What this covers

- The whole of `python/`: the `kpctoolbox` facade, the `maplib_*` example MAP
  library, the demos and their trace loaders, the `_extras` primitives, the
  pytest suite, and `python/README.md`.
- The port translates an existing, published method rather than inventing one. It
  also does not reimplement the numerics: the fitting routines and numerical
  primitives are imported from the `line-solver` package and re-exported under a
  flat namespace that mirrors the MATLAB KPC-Toolbox API.
- Part of the recent maintenance of `matlab/` (bug fixes and code cleanups) was
  also AI-assisted. The git history records which commits.

## What this does not cover

The KPC method, the MATLAB reference implementation under `matlab/`, and the
results in the papers cited in `README.md` are the work of the authors listed in
`AUTHORS.txt`. `matlab/` remains the definitive implementation: where the two
languages disagree, the MATLAB behaviour is the intended one.

## How the port was checked

During development the Python results were compared against the MATLAB toolbox
across MAP algebra and statistics, the counting process, the composition kernels,
the CTMC/DTMC solvers, the trace estimators, the example MAP library, and the
bundled traces. Agreement is to machine precision apart from a small number of
deviations, each with an identified cause, listed under "Parity with the MATLAB
toolbox" in `python/README.md`. `python/tests/` holds a self-contained pytest
suite covering the example library, PH fitting and MAP fitting.

## Caveats

Generated code carries the usual risks. Passing tests are evidence, not proof,
and the checks cover the inputs they exercise rather than every input. Treat
`python/` as you would any third-party port: validate it on your own data before
relying on it, and prefer `matlab/` when you need reference behaviour. Please
report any discrepancy to the maintainer (see `CONTACT-US.TXT`).
