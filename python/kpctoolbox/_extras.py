"""Primitives the MATLAB toolbox ships that have no home in ``line-solver``.

Everything else in :mod:`kpctoolbox` is re-exported from ``line_solver``. These
three are implemented here because ``line-solver`` does not provide them, and
their absence was the only thing keeping the Python API from covering the
MATLAB one for this part of the toolbox.
"""

import warnings

import numpy as np

__all__ = ["krons", "aph_from2moments", "map2mmpp"]


def krons(A, B):
    """Kronecker sum of ``A`` and ``B``.

    Port of ``matlab/basic/krons.m``::

        krons(A, B) = kron(A, I_B) + kron(I_A, B)
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return np.kron(A, np.eye(*B.shape)) + np.kron(np.eye(*A.shape), B)


def aph_from2moments(MEAN, SCV):
    """Acyclic PH matching a target mean and SCV, in renewal MAP ``(D0, D1)`` form.

    Port of ``matlab/aph/aph_from2moments.m``: the canonical acyclic
    construction of Bobbio, Horvath and Telek.

    Args:
        MEAN: target mean
        SCV: target squared coefficient of variation (variance / mean^2)

    Returns:
        Tuple ``(D0, D1)``, where ``D0`` is the PH subgenerator ``T`` and
        ``D1 = (-T e) alpha`` restarts the phase process on each arrival.
    """
    cv2 = float(SCV)
    lam = 1.0 / float(MEAN)
    N = max(int(np.ceil(1.0 / cv2)), 2)
    p = 1.0 / (cv2 + 1.0 + (cv2 - 1.0) / (N - 1))

    T = -lam * p * N * np.eye(N)
    for i in range(N - 1):
        T[i, i + 1] = -T[i, i]
    T[N - 1, N - 1] = -lam * N

    alpha = np.zeros(N)
    alpha[0] = p
    alpha[N - 1] = 1.0 - p

    d = -T @ np.ones(N)
    return T, np.outer(d, alpha)


def map2mmpp(D0, D1):
    """Read a MAP as an MMPP: return its generator ``Q`` and rate matrix ``LAMBDA``.

    Port of ``matlab/map/map2mmpp.m``. Warns, as MATLAB does, when ``D1`` is not
    diagonal, i.e. when the MAP is not an MMPP and ``LAMBDA`` is not a rate
    matrix.

    Args:
        D0: hidden transition matrix
        D1: visible transition matrix

    Returns:
        Tuple ``(Q, LAMBDA)`` with ``Q = D0 + D1`` and ``LAMBDA = D1``.
    """
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    if np.linalg.norm(D1 - np.diag(np.diag(D1))) > 1e-10:
        warnings.warn("The MAP is not a MMPP, LAMBDA is not diagonal", stacklevel=2)
    return D0 + D1, D1
