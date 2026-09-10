# AI Disclaimer

The Python package under `python/` (PyPI `kpc-toolbox`, `import kpctoolbox`) is a
port of the MATLAB toolbox produced with Claude Opus 5 (Anthropic), used through
Claude Code.

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
