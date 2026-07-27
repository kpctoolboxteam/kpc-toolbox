"""Example MAP library (``maplib_*``).

Ready-made Markovian Arrival Processes for demos, tests, and experimentation.
Each constructor returns a MAP as a ``(D0, D1)`` tuple of numpy arrays, matching
the representation used throughout ``line_solver`` and :mod:`kpctoolbox`.

These mirror the MATLAB KPC-Toolbox ``examples/maplib_*.m`` files. They are
provided here because ``line-solver`` does not ship this example library.
"""

from importlib import resources

import numpy as np

from kpctoolbox import (
    map_erlang,
    map_exponential,
    map_hyperexp,
    map_mmpp2,
    map_mean,
)

__all__ = [
    "maplib_exp",
    "maplib_erlang2",
    "maplib_erlang3",
    "maplib_hyper4",
    "maplib_hyper25",
    "maplib_hmap4",
    "maplib_hmap25",
    "maplib_circul4",
    "maplib_saw",
    "maplib_bcaug89",
]


def _scale_to_mean(D0, D1, new_mean):
    """Rescale a feasible MAP to a target mean inter-arrival time.

    Reproduces MATLAB ``map_scale(MAP, NEWMEAN)``: multiply both matrices by
    ``map_mean(MAP)/NEWMEAN`` (scaling the time axis) so the mean becomes
    ``NEWMEAN`` while normalized moments and correlations are preserved.
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    ratio = map_mean(D0, D1) / new_mean
    return D0 * ratio, D1 * ratio


def maplib_exp():
    """Poisson process, mean = 1."""
    return map_exponential(1.0)


def maplib_erlang2():
    """Erlang-2 process, mean = 1."""
    return map_erlang(1.0, 2)


def maplib_erlang3():
    """Erlang-3 process, mean = 1."""
    return map_erlang(1.0, 3)


def maplib_hyper4():
    """Hyper-exponential process, mean = 1, scv = 4, p = 0.99, no correlations."""
    return map_hyperexp(1.0, 4.0, 0.99)


def maplib_hyper25():
    """Hyper-exponential process, mean = 1, scv = 25, p = 0.99, no correlations."""
    return map_hyperexp(1.0, 25.0, 0.99)


def maplib_hmap4():
    """Hyper-exponential MMPP(2), mean = 1, scv = 4, weak correlations (acf1 ~ 0.1)."""
    return map_mmpp2(1.0, 4.0, -1, 0.1)


def maplib_hmap25():
    """Hyper-exponential MMPP(2), mean = 1, scv = 25, strong correlations (acf1 ~ 0.48)."""
    return map_mmpp2(1.0, 25.0, -1, -1)


def maplib_circul4():
    """MAP(4) with circulant embedded process (complex eigenvalues, oscillating ACF)."""
    D0 = np.array([
        [-1.969053238281377, 0.0, 0.0, 0.0],
        [0.0, -11.693745881208926, 0.0, 0.0],
        [0.0, 0.0, -3.809781645410352, 0.0],
        [0.0, 0.0, 0.0, -1.248416660038255],
    ])
    D1 = np.array([
        [0.057536282157225, 0.036674063864832, 0.099322454880930, 1.775520437378390],
        [10.544399916640993, 0.341694500390130, 0.217798673458358, 0.589852790719446],
        [0.192171983075854, 3.435328736606019, 0.111322877130061, 0.070958048598418],
        [0.023252043890964, 0.062972298045875, 1.125713236781891, 0.036479081319525],
    ])
    return _scale_to_mean(D0, D1, 1.0)


def maplib_saw():
    """MAP(2) with saw-like autocorrelation (-1 eigenvalue in embedded process)."""
    D0 = np.array([[-1.0, 0.0], [0.0, -2.0]])
    D1 = np.array([[0.0, 1.0], [2.0, 0.0]])
    return _scale_to_mean(D0, D1, 1.0)


def maplib_bcaug89():
    """MAP(16) obtained by KPC fitting of the Bellcore Aug89 trace.

    Reference: G. Casale, E. Zhang, E. Smirni, "Interarrival times
    characterization and fitting for Markovian traffic analysis", 2008,
    W&M Tech. Rep. WM-CS-2008-02.
    """
    with resources.files("kpctoolbox").joinpath("data/maplib_bcaug89.npz").open("rb") as fh:
        data = np.load(fh)
        return data["D0"], data["D1"]
