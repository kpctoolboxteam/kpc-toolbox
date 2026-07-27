"""Tests for the standalone kpctoolbox package.

These exercise the facade over line-solver plus the net-new example library and
demos. They require ``line-solver`` to be importable (``pip install line-solver``
or PYTHONPATH). Mirrors the MATLAB KPC-Toolbox ``tests/`` scripts, converted to
assertions.
"""

import warnings

import numpy as np
import pytest

import kpctoolbox as k
from kpctoolbox import examples as ex


def _moments(D0, D1, n=3):
    return [float(np.real(k.map_moment(D0, D1, i))) for i in range(1, n + 1)]


def _best_rel_err(PH, E):
    best = None
    for ent in PH:
        (a, b) = ent[0]
        got = _moments(a, b, len(E))
        rel = max(abs(g - t) / abs(t) for g, t in zip(got, E) if t != 0)
        best = rel if best is None else min(best, rel)
    return best


# --------------------------------------------------------------------------- #
# Example MAP library (maplib_*)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name,order,scv,acf1", [
    ("maplib_exp", 1, 1.0, 0.0),
    ("maplib_erlang2", 2, 0.5, 0.0),
    ("maplib_erlang3", 3, 1.0 / 3.0, 0.0),
    ("maplib_hyper4", 2, 4.0, 0.0),
    ("maplib_hyper25", 2, 25.0, 0.0),
    ("maplib_hmap4", 2, 4.0, 0.1),
    ("maplib_hmap25", 2, 25.0, 0.48),
])
def test_maplib_moments_and_correlation(name, order, scv, acf1):
    D0, D1 = getattr(ex, name)()
    assert D0.shape[0] == order
    assert k.map_isfeasible(D0, D1)
    assert float(k.map_mean(D0, D1)) == pytest.approx(1.0, rel=1e-6)
    assert float(k.map_scv(D0, D1)) == pytest.approx(scv, rel=1e-4)
    assert float(np.ravel(k.map_acf(D0, D1, 1))[0]) == pytest.approx(acf1, abs=1e-3)


@pytest.mark.parametrize("name,order", [
    ("maplib_circul4", 4), ("maplib_saw", 2), ("maplib_bcaug89", 16)])
def test_maplib_feasible(name, order):
    D0, D1 = getattr(ex, name)()
    assert D0.shape[0] == order
    assert k.map_isfeasible(D0, D1)


# --------------------------------------------------------------------------- #
# PH fitting (kpcfit_ph_auto) — ports of tests/test_kpcfit_ph_*.m
# --------------------------------------------------------------------------- #
def _hyper2():
    H0 = np.array([[-12.000975230358, 0.000975230358],
                   [0.000080862578, -0.088000599518]])
    H1 = np.array([[12.0, 0.0], [0.0, 0.08791973694]])
    return k.map_renewal(H0, H1) if False else (H0, H1)


def _hypo2():
    D0 = np.array([[-2.389377046831244, 1.907112341413688],
                   [0.1141504016723399, -4.480622953168757]])
    D1 = np.array([[0.0, 0.4822647054175562],
                   [4.3495615093454, 0.01691104215101627]])
    return D0, D1


@pytest.mark.parametrize("maker", [
    pytest.param(lambda: _hyper2(), id="hyper2"),
    pytest.param(lambda: _hypo2(), id="hypo2"),
    pytest.param(lambda: k.map_erlang(2.0, 2), id="erlang2"),
])
def test_kpcfit_ph_auto_matches_three_moments(maker):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D0, D1 = maker()
        E = np.array(_moments(D0, D1, 3))
        opts = k.kpcfit_ph_options(E, verbose=False)
        PH = k.kpcfit_ph_auto(E, opts)
    assert len(PH) >= 1
    assert _best_rel_err(PH, E) < 1e-6


def test_kpcfit_ph_exact_erlang3_needs_seven_moments():
    # SCV = 1/3 requires >= 3 phases; with 7 moments the exact low-variability
    # path recognizes the Erlang-3 moment set and fits it exactly. (The full
    # kpcfit_ph_auto also runs a much slower approximate KPC search at 4 states,
    # so we exercise the exact method here.)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D0, D1 = k.map_erlang(3.0, 3)
        E = np.array(_moments(D0, D1, 7))
        opts = k.kpcfit_ph_options(E, verbose=False)
        PH = k.kpcfit_ph_exact(E, opts)
    assert len(PH) >= 1
    # kpcfit_ph_exact returns a list of (D0, D1) pairs directly
    assert _best_rel_err([(pair, 0) for pair in PH], E[:3]) < 1e-6


def test_kpcfit_ph_auto_preserves_non_unit_mean():
    # regression: fitted PH must keep the target mean, not collapse to 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D0, D1 = k.map_hyperexp(5.0, 8.0, 0.6)
        E = np.array(_moments(D0, D1, 3))
        PH = k.kpcfit_ph_auto(E, k.kpcfit_ph_options(E, verbose=False))
    assert len(PH) >= 1
    (a, b) = PH[0][0]
    assert float(np.real(k.map_mean(a, b))) == pytest.approx(5.0, rel=1e-6)


# --------------------------------------------------------------------------- #
# MAP fitting (kpcfit_init -> kpcfit_auto) on a bundled trace subsample
# --------------------------------------------------------------------------- #
def test_kpcfit_auto_on_trace_subsample():
    from kpctoolbox import demo
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        S = demo.load_trace("BCAUG89")[:15000]
        trace = k.kpcfit_init(S)
        res = k.kpcfit_auto(trace, OnlyAC=True, NumStates=2, NumMAPs=1,
                            MaxRunsAC=1, MaxIterAC=20)
    MAP = res[0] if isinstance(res, tuple) else res
    D0, D1 = MAP
    assert k.map_isfeasible(D0, D1)
    assert float(k.map_mean(D0, D1)) > 0
