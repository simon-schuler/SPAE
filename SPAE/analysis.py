"""
Post-processing for SPAE MCMC runs — replaces the hand-copied analysis
notebook (see Verification/SPAE_test/spae_plots.ipynb).

Loading
-------
load_sampler        — unpickle a saved emcee sampler
load_flat_blob       — load the convenience .npy blob cache (cross-check only)

Numeric core
------------
clip_run              — burn-in + acceptance-fraction walker filtering
median_solution        — median (Teff, logg, [Fe/H], micro) + per-species table
balance_stats           — EP/REW slope diagnostics for a set of abundances
traditional_solution     — excitation/ionization-balance solution from the posterior
refine_traditional_solution — optional single-best-sample refinement (costly)
analyze_run                — one-call convenience chaining all of the above

I/O
---
summarize            — human-readable text summary of an AnalysisResult

Plotting (all return (fig, ax); matplotlib imported lazily; no rcParams set)
--------
plot_trace, plot_balance, plot_corner, plot_species_histogram
"""
from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from scipy.stats import linregress

from .utils import sigma_func
from .abunds import abunds_func, rel_abunds, abs_abunds


_SPECIES_DTYPE = np.dtype([
    ('species', 'U10'), ('n_lines', 'i4'), ('mean', 'f8'),
    ('std', 'f8'), ('sigma_mean', 'f8'), ('tot_uncert', 'f8'),
])


# --------------------------------------------------------------------------- #
# Result containers
# --------------------------------------------------------------------------- #

@dataclass
class ParamEstimate:
    """One stellar parameter's median and 1-sigma offsets, from sigma_func()."""
    median: float
    lower: float   # 16th percentile - median (<= 0)
    upper: float   # 84th percentile - median (>= 0)


@dataclass
class ClippedRun:
    flat_chain: np.ndarray        # (n_kept, 4)
    flat_blob: np.ndarray         # (n_kept,) structured, same dtype as sampler blobs
    n_walkers_kept: int
    n_walkers_total: int
    burn_in: int
    a_frac: float


@dataclass
class MedianSolution:
    teff: ParamEstimate
    logg: ParamEstimate
    feh: ParamEstimate
    micro: ParamEstimate
    sol: Tuple[float, float, float, float]
    el_found: List[str]
    abundances: List[np.ndarray]     # abunds_func() output at `sol`
    species_table: np.ndarray        # dtype=_SPECIES_DTYPE, one row per el_found


@dataclass
class BalanceStats:
    ep: object    # scipy.stats.linregress result for Fe I abund vs EP
    rew: object   # scipy.stats.linregress result for Fe I abund vs logRW
    fe1: np.ndarray
    fe2: np.ndarray


@dataclass
class TraditionalSolution:
    max_cor: float
    n_checked: int
    n_accepted: int
    n_unique: int
    uniq_params: np.ndarray   # (n_unique, 4)
    mean: Tuple[float, float, float, float]
    std: Tuple[float, float, float, float]


@dataclass
class RefinedTraditionalSolution:
    p_threshold: float
    n_qualifying: int
    qualifying_params: np.ndarray
    mean: Tuple[float, float, float, float]
    std: Tuple[float, float, float, float]
    best_params: Tuple[float, float, float, float]
    best_ep: object
    best_rew: object


@dataclass
class AnalysisResult:
    clipped: ClippedRun
    median: MedianSolution
    median_balance: Optional[BalanceStats]
    traditional: TraditionalSolution
    refined: Optional[RefinedTraditionalSolution] = None


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

def load_sampler(pickle_path: str):
    """Unpickle a saved emcee sampler (as written by pickle.dump(sampler, f))."""
    with open(pickle_path, 'rb') as f:
        return pickle.load(f)


def load_flat_blob(npy_path: str) -> np.ndarray:
    """
    Load the convenience `{star}_flat_blob_{run}.npy` cache.

    This is the UNCLIPPED sampler.get_blobs().reshape(-1) — no burn-in or
    acceptance-fraction filtering applied, and not used by any function
    below. clip_run(sampler, ...) is the correct input for all analysis
    in this module; this loader exists only for quick manual inspection.
    """
    return np.load(npy_path)


# --------------------------------------------------------------------------- #
# Numeric core
# --------------------------------------------------------------------------- #

def clip_run(sampler, burn_in: int = 100, a_frac: float = 0.40) -> ClippedRun:
    """
    Discard low-acceptance walkers and the first `burn_in` steps.

    Keeps walker i only if sampler.acceptance_fraction[i] > a_frac; for kept
    walkers, concatenates chain[burn_in:, i, :] / blobs[burn_in:, i].
    """
    chain = sampler.get_chain()
    blobs = sampler.get_blobs()
    n_steps, n_walkers, n_variables = chain.shape

    flat_chain = None
    flat_blob = None
    n_kept = 0
    for i in range(n_walkers):
        if sampler.acceptance_fraction[i] > a_frac:
            if flat_chain is None:
                flat_chain = chain[burn_in:, i, :]
                flat_blob = blobs[burn_in:, i]
            else:
                flat_chain = np.append(flat_chain, chain[burn_in:, i, :], axis=0)
                flat_blob = np.append(flat_blob, blobs[burn_in:, i])
            n_kept += 1

    if flat_chain is None:
        raise ValueError(
            f"No walkers exceeded acceptance fraction {a_frac} "
            f"(max was {np.max(sampler.acceptance_fraction):.3f}); "
            "try a lower a_frac."
        )

    return ClippedRun(flat_chain=flat_chain, flat_blob=flat_blob,
                       n_walkers_kept=n_kept, n_walkers_total=n_walkers,
                       burn_in=burn_in, a_frac=a_frac)


def median_solution(clipped: ClippedRun, linelist) -> MedianSolution:
    """
    Median (Teff, logg, [Fe/H], micro) with 1-sigma offsets, re-derived
    abundances at that point, and a per-species summary table built from
    the already-computed MCMC blob fields (not a recompute).
    """
    t = ParamEstimate(*sigma_func(clipped.flat_chain[:, 0]))
    g = ParamEstimate(*sigma_func(clipped.flat_chain[:, 1]))
    fe = ParamEstimate(*sigma_func(clipped.flat_chain[:, 2]))
    m = ParamEstimate(*sigma_func(clipped.flat_chain[:, 3]))
    sol = (t.median, g.median, fe.median, m.median)

    el_found, abundances = abunds_func(sol, linelist)

    rows = []
    for i, el in enumerate(el_found):
        std = np.std(clipped.flat_blob[el])
        sigma_mean = np.median(clipped.flat_blob[el + '_sigma_mean'])
        tot_uncert = np.sqrt(std ** 2 + sigma_mean ** 2)
        rows.append((el, len(abundances[i]), np.mean(clipped.flat_blob[el]),
                     std, sigma_mean, tot_uncert))
    species_table = np.array(rows, dtype=_SPECIES_DTYPE)

    return MedianSolution(teff=t, logg=g, feh=fe, micro=m, sol=sol,
                           el_found=el_found, abundances=abundances,
                           species_table=species_table)


def balance_stats(el_found, abundances, sun_el=None, sun_abs=None,
                   fe1_label: str = 'Fe I ', fe2_label: str = 'Fe II ') -> BalanceStats:
    """
    EP/REW slope diagnostics for one set of abundances (absolute if
    sun_el/sun_abs are None, else differential vs. the Sun — same
    branching as abunds.obj_func).
    """
    if sun_el is None or sun_abs is None:
        fe1 = abs_abunds(el_found, abundances, fe1_label)
        fe2 = abs_abunds(el_found, abundances, fe2_label)
    else:
        fe1 = rel_abunds(el_found, abundances, sun_el, sun_abs, fe1_label)
        fe2 = rel_abunds(el_found, abundances, sun_el, sun_abs, fe2_label)

    ep = linregress(fe1['EP'], fe1['abund'])
    rew = linregress(fe1['logRWin'], fe1['abund'])
    return BalanceStats(ep=ep, rew=rew, fe1=fe1, fe2=fe2)


def traditional_solution(clipped: ClippedRun, max_cor: float = 0.01,
                          fe1_label: str = 'Fe I ', fe2_label: str = 'Fe II '
                          ) -> TraditionalSolution:
    """
    The excitation/ionization-balance ("traditional") solution: posterior
    samples with |ep_r| < max_cor, |rew_r| < max_cor, and matching Fe I/Fe II
    (rounded to 2 decimals), deduplicated and averaged.
    """
    flat_blob = clipped.flat_blob
    flat_chain = clipped.flat_chain
    n_checked = len(flat_blob)

    accepted = []
    for i in range(n_checked):
        f1 = np.round(np.mean(flat_blob[i][fe1_label]), 2)
        f2 = np.round(np.mean(flat_blob[i][fe2_label]), 2)
        if (abs(flat_blob[i]['ep_r']) < max_cor
                and abs(flat_blob[i]['rew_r']) < max_cor
                and f1 == f2):
            accepted.append(flat_chain[i, :])

    n_accepted = len(accepted)
    if n_accepted == 0:
        uniq = np.empty((0, flat_chain.shape[1]))
        mean = tuple(np.full(flat_chain.shape[1], np.nan))
        std = tuple(np.full(flat_chain.shape[1], np.nan))
    else:
        arr = np.array(accepted).reshape(n_accepted, flat_chain.shape[1])
        uniq = np.unique(arr, axis=0)
        mean = tuple(uniq.mean(axis=0))
        std = tuple(uniq.std(axis=0))

    return TraditionalSolution(max_cor=max_cor, n_checked=n_checked,
                                n_accepted=n_accepted, n_unique=len(uniq),
                                uniq_params=uniq, mean=mean, std=std)


def refine_traditional_solution(clipped: ClippedRun, traditional: TraditionalSolution,
                                 linelist, sun_el=None, sun_abs=None,
                                 p_threshold: float = 0.94,
                                 fe1_label: str = 'Fe I ', fe2_label: str = 'Fe II '
                                 ) -> Optional[RefinedTraditionalSolution]:
    """
    Optional, costlier refinement: re-derive abundances (one MOOG call per
    unique traditional_solution() sample) and keep only those whose EP/REW
    slopes are consistent with zero at `p_threshold`; report their mean/std
    plus the single sample maximizing ep_p**2 + rew_p**2.

    Returns None if there are no unique samples, or none qualify.
    """
    if traditional.n_unique == 0:
        return None

    qualifying = []
    stats_list = []
    for row in traditional.uniq_params:
        el_found, abundances = abunds_func(tuple(row), linelist)
        bstats = balance_stats(el_found, abundances, sun_el, sun_abs,
                                fe1_label, fe2_label)
        if bstats.ep.pvalue > p_threshold and bstats.rew.pvalue > p_threshold:
            qualifying.append(row)
            stats_list.append(bstats)

    if not qualifying:
        return None

    qualifying = np.array(qualifying)
    scores = [b.ep.pvalue ** 2 + b.rew.pvalue ** 2 for b in stats_list]
    best_idx = int(np.argmax(scores))

    return RefinedTraditionalSolution(
        p_threshold=p_threshold, n_qualifying=len(qualifying),
        qualifying_params=qualifying,
        mean=tuple(qualifying.mean(axis=0)), std=tuple(qualifying.std(axis=0)),
        best_params=tuple(qualifying[best_idx]),
        best_ep=stats_list[best_idx].ep, best_rew=stats_list[best_idx].rew,
    )


def analyze_run(sampler, linelist, sun_el=None, sun_abs=None,
                 burn_in: int = 100, a_frac: float = 0.40, max_cor: float = 0.01,
                 refine_p_threshold: Optional[float] = None) -> AnalysisResult:
    """
    One-call convenience: clip_run -> median_solution -> balance_stats(at
    median) -> traditional_solution -> (refine_traditional_solution only if
    refine_p_threshold is given, since it costs one MOOG call per unique
    balance sample).
    """
    clipped = clip_run(sampler, burn_in=burn_in, a_frac=a_frac)
    median = median_solution(clipped, linelist)
    median_balance = balance_stats(median.el_found, median.abundances, sun_el, sun_abs)
    traditional = traditional_solution(clipped, max_cor=max_cor)

    refined = None
    if refine_p_threshold is not None:
        refined = refine_traditional_solution(
            clipped, traditional, linelist, sun_el, sun_abs,
            p_threshold=refine_p_threshold)

    return AnalysisResult(clipped=clipped, median=median,
                           median_balance=median_balance,
                           traditional=traditional, refined=refined)


def summarize(result: AnalysisResult) -> str:
    """
    Human-readable text summary of an AnalysisResult: clipping info, MCMC
    median solution, EP/REW slopes at the median, per-species [X/H] table,
    and the traditional (excitation/ionization-balance) solution.

    The traditional-solution block reports the MEAN of all unique posterior
    samples satisfying the balance criteria — not a single "best" sample.
    If `result.refined` is set (i.e. analyze_run() was called with
    refine_p_threshold), a second block reports that refinement's mean of
    only the p-qualifying samples, plus the single best sample among them
    (the one maximizing ep_p**2 + rew_p**2) — clearly labeled as distinct
    from the plain traditional-solution mean above it.
    """
    c, m, b, t = result.clipped, result.median, result.median_balance, result.traditional
    lines = []

    lines.append(f"Clipping: burn_in={c.burn_in}, acceptance_fraction > {c.a_frac} "
                 f"-> {c.n_walkers_kept}/{c.n_walkers_total} walkers kept "
                 f"({len(c.flat_chain)} posterior samples)")
    lines.append("")

    lines.append("=== MCMC median solution ===")
    lines.append(f"Teff  = {m.teff.median:.1f} ({m.teff.lower:+.1f}/{m.teff.upper:+.1f}) K")
    lines.append(f"logg  = {m.logg.median:.3f} ({m.logg.lower:+.3f}/{m.logg.upper:+.3f})")
    lines.append(f"[Fe/H]= {m.feh.median:.3f} ({m.feh.lower:+.3f}/{m.feh.upper:+.3f})")
    lines.append(f"micro = {m.micro.median:.3f} ({m.micro.lower:+.3f}/{m.micro.upper:+.3f}) km/s")
    lines.append("")

    if b is not None:
        lines.append("=== Median-solution EP/REW slopes ===")
        lines.append(f"EP  slope = {b.ep.slope:+.5f} (r={b.ep.rvalue:+.4f}, p={b.ep.pvalue:.4f})")
        lines.append(f"REW slope = {b.rew.slope:+.5f} (r={b.rew.rvalue:+.4f}, p={b.rew.pvalue:.4f})")
        lines.append("")

    lines.append("=== [X/H] by species ===")
    for row in m.species_table:
        lines.append(f"{row['species'].strip():6s} n={row['n_lines']:3d}  "
                     f"[X/H]={row['mean']:+.3f}  sigma={row['std']:.3f}  "
                     f"tot_uncert={row['tot_uncert']:.3f}")
    lines.append("")

    lines.append("=== Traditional (excitation/ionization balance) solution ===")
    lines.append(f"MEAN of {t.n_unique} unique posterior samples satisfying "
                 f"|EP_r|,|REW_r| < {t.max_cor} and Fe I = Fe II "
                 f"(out of {t.n_checked} samples checked) -- not a single best sample")
    lines.append(f"Teff  = {t.mean[0]:.1f} +/- {t.std[0]:.1f} K")
    lines.append(f"logg  = {t.mean[1]:.3f} +/- {t.std[1]:.3f}")
    lines.append(f"[Fe/H]= {t.mean[2]:.3f} +/- {t.std[2]:.3f}")
    lines.append(f"micro = {t.mean[3]:.3f} +/- {t.std[3]:.3f} km/s")

    if result.refined is not None:
        r = result.refined
        lines.append("")
        lines.append(f"=== Refined traditional solution (p > {r.p_threshold}) ===")
        lines.append(f"MEAN of {r.n_qualifying} qualifying samples:")
        lines.append(f"Teff  = {r.mean[0]:.1f} +/- {r.std[0]:.1f} K")
        lines.append(f"logg  = {r.mean[1]:.3f} +/- {r.std[1]:.3f}")
        lines.append(f"[Fe/H]= {r.mean[2]:.3f} +/- {r.std[2]:.3f}")
        lines.append(f"micro = {r.mean[3]:.3f} +/- {r.std[3]:.3f} km/s")
        lines.append(
            f"SINGLE BEST sample (max ep_p^2+rew_p^2): "
            f"Teff={r.best_params[0]:.1f} K, logg={r.best_params[1]:.3f}, "
            f"[Fe/H]={r.best_params[2]:.3f}, micro={r.best_params[3]:.3f} km/s"
        )

    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Plotting — (fig, ax) return, matplotlib imported lazily, no rcParams set
# --------------------------------------------------------------------------- #

_DEFAULT_PARAM_LABELS = ['Teff (K)', 'log g', '[Fe/H]', 'micro (km/s)']


def plot_trace(chain: np.ndarray, labels=None, figsize=(10, 7), axes=None):
    """
    Trace plot: one stacked subplot per parameter, chain value vs step.

    `chain` is the FULL unclipped sampler.get_chain(), shape
    (n_steps, n_walkers, 4) — burn-in is a visual judgment call made by
    looking at this plot, so it intentionally isn't pre-clipped.
    """
    import matplotlib.pyplot as plt

    n_steps, n_walkers, n_dim = chain.shape
    if labels is None:
        labels = _DEFAULT_PARAM_LABELS
    if axes is None:
        fig, axes = plt.subplots(n_dim, figsize=figsize, sharex=True)
    else:
        fig = axes[0].get_figure()

    for i in range(n_dim):
        axes[i].plot(chain[:, :, i], 'k', alpha=0.3)
        axes[i].set_xlim(0, n_steps)
        axes[i].set_ylabel(labels[i], fontsize=12)
    axes[-1].set_xlabel('Step Number', fontsize=12)

    fig.tight_layout()
    return fig, axes


def plot_balance(balance: BalanceStats, ylims=None, figsize=(12, 9), axes=None):
    """EP and REW correlation scatter + linear fit, 2 stacked subplots."""
    import matplotlib.pyplot as plt

    if axes is None:
        fig, axes = plt.subplots(2, 1, figsize=figsize)
    else:
        fig = axes[0].get_figure()

    fe1 = balance.fe1
    axes[0].scatter(fe1['EP'], fe1['abund'], marker='.')
    axes[0].plot(fe1['EP'], balance.ep.intercept + balance.ep.slope * fe1['EP'],
                 color='black')
    axes[0].set_xlabel('EP (eV)', fontsize=12)
    axes[0].set_ylabel('Fe I abundance', fontsize=12)

    axes[1].scatter(fe1['logRWin'], fe1['abund'], marker='.')
    axes[1].plot(fe1['logRWin'],
                 balance.rew.intercept + balance.rew.slope * fe1['logRWin'],
                 color='black')
    axes[1].set_xlabel('log RW', fontsize=12)
    axes[1].set_ylabel('Fe I abundance', fontsize=12)

    if ylims is not None:
        axes[0].set_ylim(ylims)
        axes[1].set_ylim(ylims)

    fig.tight_layout()
    return fig, axes


def plot_corner(clipped: ClippedRun, truths=None, ranges=None, labels=None,
                 quantiles=(0.16, 0.84), levels=None, bins=30, fig=None):
    """
    Corner plot of the clipped posterior. Requires the optional `corner`
    package (pip install corner) — imported lazily, only when called.
    """
    try:
        import corner
    except ImportError:
        raise ImportError(
            "corner is required for plot_corner. Install it with:  pip install corner"
        )

    if labels is None:
        labels = _DEFAULT_PARAM_LABELS
    if levels is None:
        levels = (1 - np.exp(-0.5),)

    fig = corner.corner(clipped.flat_chain, truths=truths, range=ranges,
                         labels=labels, levels=levels,
                         quantiles=list(quantiles), bins=bins, fig=fig)
    return fig, np.array(fig.axes)


def plot_species_histogram(flat_blob: np.ndarray, species_a: str, species_b: str,
                            labels=None, colors=('orange', 'dodgerblue'),
                            density: bool = True, ax=None):
    """
    Overlaid histograms of two species' per-sample abundance blob fields
    (e.g. species_a='Fe I ', species_b='Fe II ' — note MOOG-style trailing
    space in the label, matching the blob field names).
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    if labels is None:
        labels = (species_a.strip(), species_b.strip())

    ax.hist(flat_blob[species_a], histtype='step', color=colors[0],
            density=density, label=labels[0])
    ax.hist(flat_blob[species_b], histtype='step', color=colors[1],
            density=density, label=labels[1])
    ax.set_xlabel('log N', fontsize=12)
    ax.set_ylabel('(Normalized) Number of Lines' if density else 'Number of Lines',
                  fontsize=12)
    ax.legend(fontsize=10)

    fig.tight_layout()
    return fig, ax
