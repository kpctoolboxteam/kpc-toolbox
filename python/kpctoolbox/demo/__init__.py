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

__all__ = [
    "load_trace",
    "demo_kpcfit_bcaug89",
    "demo_kpcfit_bcaug89_thirdord",
    "demo_kpcfit_dec_pkt",
    "demo_kpcfit_ph",
    "demo_run",
]

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


def demo_kpcfit_ph(trace_name="DEC-PKT-1-UDP", n_moments=11, states=4, runs=3,
                   plot=False):
    """Fit PH distributions to the first moments of a bundled trace.

    Mirrors ``demo/demo_kpcfit_ph_dec_pkt_1_udp.m`` and
    ``demo/demo_kpcfit_ph_livemaps.m``, which go through
    ``tests/test_kpcfit_ph_run.m``: eleven raw moments and a fixed PH order.

    ``n_moments`` has to be at least ``2*states - 1`` or the fitter reduces the
    order to what the moments can support, which for a trace with SCV < 1/2 can
    leave it with nothing to return.

    Returns the list of fitted PH results from :func:`kpctoolbox.kpcfit_ph_auto`.
    """
    S = load_trace(trace_name)
    E = np.array([np.mean(S ** i) for i in range(1, n_moments + 1)])
    options = k.kpcfit_ph_options(E, verbose=False, runs=runs,
                                  min_num_states=states, max_num_states=states)
    PH = k.kpcfit_ph_auto(E, options)

    print(f"{trace_name}: {len(S)} samples, PH order {states}, "
          f"target E[X^1..{n_moments}] = "
          f"{[float(x) for x in np.round(E[:3], 6)]} ...")
    print(f"  kpcfit_ph_auto returned {len(PH)} PH distribution(s)")
    for j, ent in enumerate(PH):
        (a, b), dist, _, method = ent
        fit = [float(np.real(k.map_moment(a, b, i))) for i in range(1, 4)]
        # A candidate whose mean collapses to zero or a non-finite value is a
        # failed fit that the search kept; MATLAB's kpcfit_ph_summary ranks by
        # distance and would never pick it. Say so rather than print zeros.
        ok = np.isfinite(fit[0]) and fit[0] > 0
        flag = "" if ok else "   <-- degenerate, ignore"
        print(f"    PH[{j}] order={a.shape[0]} dist={float(dist):.6g} "
              f"{method}{flag}")
        print(f"           fitE[1..3]={[float(x) for x in np.round(fit, 6)]}")

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


def demo_kpcfit_bcaug89_thirdord(max_runs_ac=10, max_iter_bc=5, max_runs_bc=1,
                                 plot=False):
    """Fit BCAUG89 to second AND third order, and compare the bicovariances.

    Mirrors ``demo/demo_kpcfit_bcaug89_thirdord.m``: an autocorrelation-only fit
    first, then a full fit that also matches the bicovariance surface, so the
    two can be held side by side. Returns ``(MAP_ord2, MAP_ord3)``.
    """
    S = load_trace("BCAUG89")
    trace = k.kpcfit_init(S)

    ord2 = k.kpcfit_auto(trace, OnlyAC=True, MaxIterAC=50)
    MAP2 = ord2[0] if isinstance(ord2, tuple) else ord2
    ord3 = k.kpcfit_auto(trace, MaxRunsAC=max_runs_ac,
                         MaxIterBC=max_iter_bc, MaxRunsBC=max_runs_bc)
    MAP3 = ord3[0] if isinstance(ord3, tuple) else ord3

    print(f"BCAUG89 third-order fit: {len(S)} samples")
    for label, M in (("second-order (AC only)", MAP2), ("third-order (AC+BC)", MAP3)):
        D0, D1 = M[0], M[1]
        bc = np.array([float(np.real(k.map_joint(D0, D1, lag, [1, 1, 1])))
                       for lag in np.asarray(trace.BCLags)])
        err = float(np.max(np.abs(bc - np.ravel(trace.BC))))
        print(f"  {label}: order {D0.shape[0]}, max |BC - BC_fit| = {err:.6e}")

    if plot:
        import matplotlib.pyplot as plt
        g = np.log10(np.ravel(trace.BCGridLags))
        X, Y = np.meshgrid(g, g)
        n = len(g)
        fig = plt.figure(figsize=(11, 4))
        for i, (label, M) in enumerate((("second-order", MAP2),
                                        ("third-order", MAP3))):
            D0, D1 = M[0], M[1]
            est = np.array([float(np.real(k.map_joint(D0, D1, lag, [1, 1, 1])))
                            for lag in np.asarray(trace.BCLags)]).reshape(n, n)
            ax = fig.add_subplot(1, 2, i + 1, projection="3d")
            ax.plot_surface(X, Y, np.ravel(trace.BC).reshape(n, n), alpha=0.5)
            ax.plot_surface(X, Y, est, alpha=0.5)
            ax.set_title(f"{label} fit")
        fig.tight_layout()
        plt.show()

    return MAP2, MAP3


def demo_kpcfit_dec_pkt(num_states=8, smooth=2, max_iter_ac=50, plot=False):
    """Fit a MAP of a given order to the DEC-PKT-1-UDP trace.

    Mirrors ``demo/demo_kpcfit_dec_pkt_1_udp8.m`` and its 16-state sibling:
    pass ``num_states=8`` or ``num_states=16``. Returns the fitted ``(D0, D1)``.
    """
    S = load_trace("DEC-PKT-1-UDP")
    trace = k.kpcfit_init(S, smooth=smooth)  # line spells this one snake_case
    result = k.kpcfit_auto(trace, OnlyAC=True, MaxRunsAC=1,
                           MaxIterAC=max_iter_ac, NumStates=num_states)
    MAP = result[0] if isinstance(result, tuple) else result
    D0, D1 = MAP[0], MAP[1]

    print(f"DEC-PKT-1-UDP: {len(S)} samples, fitted MAP order {D0.shape[0]} "
          f"(requested {num_states})")
    print(f"  trace mean={k.trace_mean(S):.6g} scv={k.trace_scv(S):.4f}")
    print(f"  MAP   mean={float(k.map_mean(D0, D1)):.6g} "
          f"scv={float(k.map_scv(D0, D1)):.4f}")

    if plot:
        import matplotlib.pyplot as plt
        lags = np.arange(1, 101)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(lags, k.trace_acf(S, lags), "k-", label="trace")
        ax[0].plot(trace.ACLags, trace.AC, "r*", label="fitted lags")
        ax[0].plot(lags, np.ravel(k.map_acf(D0, D1, lags)), "b-", label="MAP")
        ax[0].set_xlabel("Lag k"); ax[0].set_ylabel(r"ACF $\rho_k$"); ax[0].legend()
        ax[1].loglog(trace.ACLags, trace.AC, "r-")
        ax[1].loglog(trace.ACLags, np.ravel(k.map_acf(D0, D1, trace.ACLags)), "b-")
        ax[1].set_xlabel("Lag k [log]"); ax[1].set_ylabel(r"ACF $\rho_k$ [log]")
        fig.suptitle(f"DEC-PKT-1-UDP, {num_states}-state MAP fit")
        fig.tight_layout()
        plt.show()

    return D0, D1


def demo_run(plot=False):
    """Run every demo in sequence, as ``demo/demo_run.m`` does.

    MATLAB pauses between scripts; here each demo simply runs to completion and
    its result is collected. Returns a dict keyed by demo name.
    """
    out = {}
    for name, fn in (
        ("demo_kpcfit_bcaug89", lambda: demo_kpcfit_bcaug89(plot=plot)),
        ("demo_kpcfit_bcaug89_thirdord", lambda: demo_kpcfit_bcaug89_thirdord(plot=plot)),
        ("demo_kpcfit_dec_pkt_8", lambda: demo_kpcfit_dec_pkt(8, plot=plot)),
        ("demo_kpcfit_dec_pkt_16", lambda: demo_kpcfit_dec_pkt(16, plot=plot)),
        ("demo_kpcfit_ph_dec_pkt", lambda: demo_kpcfit_ph("DEC-PKT-1-UDP", plot=plot)),
        ("demo_kpcfit_ph_livemaps", lambda: demo_kpcfit_ph("LIVEMAPS-1030", states=4,
                                                           runs=3, plot=plot)),
    ):
        print("\n" + "=" * 70 + f"\n{name}\n" + "=" * 70)
        out[name] = fn()
    return out
