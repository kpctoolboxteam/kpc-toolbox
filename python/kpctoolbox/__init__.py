"""KPC-Toolbox (Python edition).

Fitting of Markovian Arrival Processes (MAPs) and phase-type (PH)
distributions to empirical traces via Kronecker Product Composition (KPC).

This package is a thin standalone distribution: it does NOT reimplement any
algorithm that already ships in ``line-solver``. Every fitting routine and
numerical primitive is imported from ``line_solver`` and re-exported here under
a flat namespace that mirrors the original MATLAB KPC-Toolbox API (``kpcfit_*``,
``map_*``, ``trace_*``, ``aph_*``, ``mmpp2_*``, ``ctmc_*``/``dtmc_*`` ...). What
this package adds on top of ``line-solver`` is the material the MATLAB toolbox
ships that has no Python home elsewhere: the example MAP library
(:mod:`kpctoolbox.examples`), the runnable demos (:mod:`kpctoolbox.demo`) with
their bundled traces, and the test suite.

Requires ``line-solver`` (>= 3.0.7): ``pip install line-solver``.

Entry points
------------
MAP fitting          : :func:`kpcfit_init` -> :func:`kpcfit_auto`
PH fitting           : :func:`kpcfit_ph_options` -> :func:`kpcfit_ph_auto`

References
----------
[1] G. Casale, E. Z. Zhang, E. Smirni. KPC-Toolbox: Best Recipes for Automatic
    Trace Fitting Using Markovian Arrival Processes. Performance Evaluation,
    67(9):873-896, 2010.
[2] G. Casale, E. Z. Zhang, E. Smirni. Trace Data Characterization and Fitting
    for Markov Modeling. Performance Evaluation, 67(2):61-79, 2010.
"""

try:
    import line_solver as _line_solver  # noqa: F401
except ImportError as _exc:  # pragma: no cover - dependency guard
    raise ImportError(
        "kpctoolbox requires the 'line-solver' package for its numerical core. "
        "Install it with: pip install line-solver"
    ) from _exc

__version__ = "0.5.0"


def version():
    """Return the KPC-Toolbox version string (parity with MATLAB kpcfit_version)."""
    return __version__


# ---------------------------------------------------------------------------
# KPC fitting engine (from line_solver.lib.kpctoolbox.kpcfit).
# Several of these are implemented in line but omitted from its public __all__;
# we import them directly from the submodule so the full MATLAB entry-point set
# is available here.
# ---------------------------------------------------------------------------
from line_solver.lib.kpctoolbox.kpcfit import (  # noqa: E402
    # Data / options / result containers
    KpcfitTraceData,
    KpcfitPhOptions,
    KpcfitResult,
    KPCFIT_TOL,
    # MAP fitting pipeline
    kpcfit_tol,
    kpcfit_init,
    kpcfit_auto,
    kpcfit_manual,
    kpcfit_sub_acfit,
    kpcfit_sub_bcfit,
    kpcfit_sub_bic,
    kpcfit_sub_compose,
    kpcfit_sub_eval_acfit,
    kpcfit_hyper_charpoly,
    # PH fitting pipeline
    kpcfit_ph_options,
    kpcfit_ph_auto,
    kpcfit_ph_exact,
    kpcfit_ph_search,
    kpcfit_ph_manual,
    kpcfit_ph_prony,
    kpcfit_ph_summary,
    logspacei,
)

# ---------------------------------------------------------------------------
# MAP algebra & statistics primitives (from line_solver.api.mam).
# ---------------------------------------------------------------------------
from line_solver.api.mam import (  # noqa: E402
    # moments / statistics
    map_moment,
    map_mean,
    map_var,
    map_scv,
    map_skew,
    map_kurt,
    map_lambda,
    map_acf,
    map_acfc,
    map_idc,
    map_gamma,
    map_gamma2,
    map_joint,
    # counting process
    map_count_mean,
    map_count_var,
    map_count_moment,
    map_count_idc,
    map_varcount,
    # constructors
    map_exponential,
    map_erlang,
    map_hyperexp,
    map_mmpp2,
    map_rand,
    map_randn,
    map_renewal,
    # composition / feasibility
    map_kpc,
    map_block,
    map_feasblock,
    map_feastol,
    map_super,
    map_sum,
    map_sumind,
    map_mixture,
    map_max,
    map_mark,
    map_stochcomp,
    map_isfeasible,
    map_normalize,
    map_scale,
    map_largemap,
    map_issym,
    # generator / stationary / transforms
    map_infgen,
    map_pie,
    map_prob,
    map_piq,
    map_embedded,
    map_timereverse,
    map_sample,
    map_cdf,
    map_pdf,
    map_pntiter,
    map_pntquad,
    # MAP(2) inverse fit
    map2_fit,
)

# map2ph lives in the lib.kpctoolbox.map module (not the api.mam layer).
from line_solver.lib.kpctoolbox.map import map2ph  # noqa: E402

# ---------------------------------------------------------------------------
# Trace statistics, APH, MMPP, Markov-chain solvers, basic utilities.
# These are re-exported wholesale from the line kpctoolbox facade.
# ---------------------------------------------------------------------------
from line_solver.lib.kpctoolbox import (  # noqa: E402
    # trace statistics
    trace_mean,
    trace_var,
    trace_scv,
    trace_acf,
    trace_skew,
    trace_joint,
    trace_bicov,
    trace_idi,
    trace_idc,
    trace_gamma,
    trace_shuffle,
    trace_iat2counts,
    trace_iat2bins,
    trace_pmf,
    trace_summary,
    autocov,
    # acyclic PH
    aph_fit,
    aph_simplify,
    aph_convpara,
    aph_convseq,
    aph_rand,
    ph2hyper,
    hyper_rand,
    ConvolutionPattern,
    # MMPP(2)
    mmpp2_fit,
    mmpp2_fit1,
    mmpp2_fit2,
    mmpp2_fit3,
    mmpp2_fit4,
    mmpp2_fitc,
    mmpp2_fitc_approx,
    mmpp2_fitc_theoretical,
    mmpp_rand,
    # CTMC / DTMC
    ctmc_solve,
    ctmc_makeinfgen,
    ctmc_timereverse,
    ctmc_uniformization,
    ctmc_transient,
    ctmc_rand,
    dtmc_solve,
    dtmc_makestochastic,
    dtmc_isfeasible,
    dtmc_rand,
    dtmc_timereverse,
    weaklyconncomp,
    # basic utilities
    minpos,
    maxpos,
    spectd,
)

__all__ = [
    "version",
    "__version__",
    # containers
    "KpcfitTraceData",
    "KpcfitPhOptions",
    "KpcfitResult",
    "KPCFIT_TOL",
    # MAP fitting pipeline
    "kpcfit_tol",
    "kpcfit_init",
    "kpcfit_auto",
    "kpcfit_manual",
    "kpcfit_sub_acfit",
    "kpcfit_sub_bcfit",
    "kpcfit_sub_bic",
    "kpcfit_sub_compose",
    "kpcfit_sub_eval_acfit",
    "kpcfit_hyper_charpoly",
    # PH fitting pipeline
    "kpcfit_ph_options",
    "kpcfit_ph_auto",
    "kpcfit_ph_exact",
    "kpcfit_ph_search",
    "kpcfit_ph_manual",
    "kpcfit_ph_prony",
    "kpcfit_ph_summary",
    "logspacei",
    # MAP primitives
    "map_moment",
    "map_mean",
    "map_var",
    "map_scv",
    "map_skew",
    "map_kurt",
    "map_lambda",
    "map_acf",
    "map_acfc",
    "map_idc",
    "map_gamma",
    "map_gamma2",
    "map_joint",
    "map_count_mean",
    "map_count_var",
    "map_count_moment",
    "map_count_idc",
    "map_varcount",
    "map_exponential",
    "map_erlang",
    "map_hyperexp",
    "map_mmpp2",
    "map_rand",
    "map_randn",
    "map_renewal",
    "map_kpc",
    "map_block",
    "map_feasblock",
    "map_feastol",
    "map_super",
    "map_sum",
    "map_sumind",
    "map_mixture",
    "map_max",
    "map_mark",
    "map_stochcomp",
    "map_isfeasible",
    "map_normalize",
    "map_scale",
    "map_largemap",
    "map_issym",
    "map_infgen",
    "map_pie",
    "map_prob",
    "map_piq",
    "map_embedded",
    "map_timereverse",
    "map_sample",
    "map_cdf",
    "map_pdf",
    "map_pntiter",
    "map_pntquad",
    "map2_fit",
    "map2ph",
    # trace statistics
    "trace_mean",
    "trace_var",
    "trace_scv",
    "trace_acf",
    "trace_skew",
    "trace_joint",
    "trace_bicov",
    "trace_idi",
    "trace_idc",
    "trace_gamma",
    "trace_shuffle",
    "trace_iat2counts",
    "trace_iat2bins",
    "trace_pmf",
    "trace_summary",
    "autocov",
    # APH
    "aph_fit",
    "aph_simplify",
    "aph_convpara",
    "aph_convseq",
    "aph_rand",
    "ph2hyper",
    "hyper_rand",
    "ConvolutionPattern",
    # MMPP(2)
    "mmpp2_fit",
    "mmpp2_fit1",
    "mmpp2_fit2",
    "mmpp2_fit3",
    "mmpp2_fit4",
    "mmpp2_fitc",
    "mmpp2_fitc_approx",
    "mmpp2_fitc_theoretical",
    "mmpp_rand",
    # CTMC / DTMC
    "ctmc_solve",
    "ctmc_makeinfgen",
    "ctmc_timereverse",
    "ctmc_uniformization",
    "ctmc_transient",
    "ctmc_rand",
    "dtmc_solve",
    "dtmc_makestochastic",
    "dtmc_isfeasible",
    "dtmc_rand",
    "dtmc_timereverse",
    "weaklyconncomp",
    # basic
    "minpos",
    "maxpos",
    "spectd",
]
