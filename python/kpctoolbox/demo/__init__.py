"""Runnable demonstrations of KPC-Toolbox fitting, with bundled traces.

Traces (interarrival-time vectors) are shipped as compressed ``.npz`` under
``kpctoolbox/data`` and loaded with :func:`load_trace`. The demo functions run
the fitting pipeline end-to-end and return the fitted model(s); pass
``plot=True`` to also render a comparison with matplotlib (optional dependency).

These mirror the MATLAB KPC-Toolbox ``demo/`` scripts; ``line-solver`` does not
ship them.
"""

from importlib import resources

import numpy as np

import kpctoolbox as k

__all__ = ["load_trace", "demo_kpcfit_bcaug89", "demo_kpcfit_ph"]

_TRACES = {
    "BCAUG89": "BCAUG89.npz",
    "DEC-PKT-1-UDP": "DEC-PKT-1-UDP.npz",
    "LIVEMAPS-1030": "LIVEMAPS-1030.npz",
}


def load_trace(name):
    """Return the interarrival-time sample vector ``S`` of a bundled trace.

    Parameters
    ----------
    name : {"BCAUG89", "DEC-PKT-1-UDP", "LIVEMAPS-1030"}
    """
    if name not in _TRACES:
        raise ValueError(f"unknown trace {name!r}; choose from {sorted(_TRACES)}")
    with resources.files("kpctoolbox").joinpath("data/" + _TRACES[name]).open("rb") as fh:
        data = np.load(fh)
        return np.asarray(data["S"], dtype=np.float64).ravel()


def demo_kpcfit_bcaug89(only_ac=True, max_iter_ac=50, plot=False):
    """Fit a MAP to the Bellcore Aug89 trace and compare autocorrelations.

    Mirrors ``demo/demo_kpcfit_bcaug89.m``. Returns the fitted MAP ``(D0, D1)``.
    """
    S = load_trace("BCAUG89")
    trace = k.kpcfit_init(S)
    result = k.kpcfit_auto(trace, OnlyAC=only_ac, MaxIterAC=max_iter_ac)
    MAP = result[0] if isinstance(result, tuple) else result
    D0, D1 = MAP[0], MAP[1]

    print(f"BCAUG89: {len(S)} samples, fitted MAP order {D0.shape[0]}")
    print(f"  trace mean={k.trace_mean(S):.6g} scv={k.trace_scv(S):.4f}")
    print(f"  MAP   mean={float(k.map_mean(D0, D1)):.6g} "
          f"scv={float(k.map_scv(D0, D1)):.4f}")

    if plot:
        import matplotlib.pyplot as plt
        lags = np.arange(1, 51)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(lags, k.trace_acf(S, lags), "k-", label="trace")
        ax[0].plot(lags, np.ravel(k.map_acf(D0, D1, lags)), "b-", label="MAP")
        ax[0].set_xlabel("Lag k"); ax[0].set_ylabel(r"ACF $\rho_k$"); ax[0].legend()
        ax[1].loglog(lags, k.trace_acf(S, lags), "k-")
        ax[1].loglog(lags, np.ravel(k.map_acf(D0, D1, lags)), "b-")
        ax[1].set_xlabel("Lag k [log]"); ax[1].set_ylabel(r"ACF $\rho_k$ [log]")
        fig.suptitle("BCAUG89 MAP fit"); fig.tight_layout()
        plt.show()

    return D0, D1


def demo_kpcfit_ph(trace_name="DEC-PKT-1-UDP", n_moments=3, plot=False):
    """Fit PH distributions to the first moments of a bundled trace.

    Mirrors ``demo/demo_kpcfit_ph_*.m``. Returns the list of fitted PH results
    from :func:`kpctoolbox.kpcfit_ph_auto`.
    """
    S = load_trace(trace_name)
    E = np.array([k.trace_moment(S, i) if hasattr(k, "trace_moment")
                  else np.mean(S ** i) for i in range(1, n_moments + 1)])
    options = k.kpcfit_ph_options(E, verbose=False)
    PH = k.kpcfit_ph_auto(E, options)

    print(f"{trace_name}: {len(S)} samples, target E[X^1..{n_moments}]="
          f"{[float(x) for x in np.round(E, 6)]}")
    print(f"  kpcfit_ph_auto returned {len(PH)} PH distribution(s)")
    for j, ent in enumerate(PH):
        (a, b) = ent[0]
        fit = [float(np.real(k.map_moment(a, b, i))) for i in range(1, n_moments + 1)]
        print(f"    PH[{j}] order={a.shape[0]} fitE={list(np.round(fit, 6))}")

    if plot and PH:
        import matplotlib.pyplot as plt
        xs = np.logspace(np.log10(max(S.min(), 1e-9)), np.log10(S.max()), 200)
        plt.figure()
        Ssort = np.sort(S)
        ccdf = 1.0 - np.arange(1, len(Ssort) + 1) / len(Ssort)
        idx = np.linspace(0, len(Ssort) - 1, min(1000, len(Ssort))).astype(int)
        plt.loglog(Ssort[idx], ccdf[idx], "k-", label="trace")
        (a, b) = PH[0][0]
        plt.loglog(xs, 1.0 - np.ravel(k.map_cdf(a, b, xs)), "-", label="best PH")
        plt.xlabel("Inter-arrival time t"); plt.ylabel("CCDF Pr(X>t)")
        plt.legend(); plt.title(trace_name); plt.show()

    return PH
