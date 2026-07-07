"""Generate the paper's figures into paper/figures/ (PDF, vector).

Reads the census (data/zeros_t5000.txt) and the density curve
(data/density_curve.txt); the phase diagram is self-contained.  Run after
`python -m primeleaf census --t-max 5000 ... --out data/zeros_t5000.txt`.

    python scripts/make_figures.py

Requires matplotlib, numpy, mpmath (and the primeleaf package on the path).
"""

from __future__ import annotations

import math
import pathlib

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mpmath as mp  # noqa: E402

from primeleaf.primes import primes_up_to  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIGS = ROOT / "paper" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

SIGMA_STAR = 2.5973851271346716
PERIOD = 2 * math.pi / math.log(2)

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.linewidth": 0.6,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
    }
)
BLUE = "#2b4b8c"
RED = "#b0303a"


def load_census(path: pathlib.Path):
    re_parts, im_parts = [], []
    for line in path.read_text().splitlines():
        a, b = line.split()
        re_parts.append(float(a))
        im_parts.append(float(b))
    order = np.argsort(im_parts)
    return np.array(re_parts)[order], np.array(im_parts)[order]


def fig_zeros(re_parts, im_parts):
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.scatter(re_parts, im_parts, s=7, color=BLUE, alpha=0.8, edgecolors="none")
    ax.axvline(SIGMA_STAR, color=RED, lw=0.8, ls="--",
               label=r"$\sigma^{*}=2.5974$")
    ax.axvline(1.0, color="gray", lw=0.8, ls=":", label=r"$\mathrm{Re}\,s=1$")
    ax.set_xlabel(r"$\mathrm{Re}\,s$")
    ax.set_ylabel(r"$\mathrm{Im}\,s$")
    ax.set_xlim(0.98, 2.62)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=8)
    fig.savefig(FIGS / "zeros.pdf")
    plt.close(fig)


def fig_staircase(im_parts, density_at_1):
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    t = np.sort(im_parts)
    counts = np.arange(1, len(t) + 1)
    ax.step(t, counts, where="post", color=BLUE, lw=1.0, label="census $N(T)$")
    tmax = t[-1]
    grid = np.linspace(0, tmax, 200)
    ax.plot(grid, density_at_1 * grid, color=RED, lw=1.0,
            label=rf"$D(1)\,T$, $D(1)={density_at_1:.4f}$")
    ax.plot(grid, grid / PERIOD, color="gray", lw=0.9, ls=":",
            label=r"one per $2\pi/\log 2$")
    ax.set_xlabel(r"$T$")
    ax.set_ylabel(r"number of zeros up to height $T$")
    ax.set_xlim(0, tmax)
    ax.set_ylim(0, len(t) * 1.02)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=8)
    fig.savefig(FIGS / "staircase.pdf")
    plt.close(fig)


def fig_spacings(im_parts):
    t = np.sort(im_parts)
    gaps = np.diff(t)
    s = gaps / gaps.mean()  # unfolded to mean 1
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.hist(s, bins=np.linspace(0, s.max() + 0.2, 26), density=True,
            color=BLUE, alpha=0.7, edgecolor="white", linewidth=0.4,
            label="tree-zero spacings")
    x = np.linspace(0, max(3.0, s.max()), 300)
    ax.plot(x, np.exp(-x), color="gray", lw=1.0, ls=":", label="Poisson $e^{-s}$")
    ax.plot(x, (np.pi / 2) * x * np.exp(-np.pi / 4 * x**2), color=RED, lw=1.0,
            label=r"GUE (Wigner)")
    ax.set_xlabel(r"normalized spacing $s$")
    ax.set_ylabel("density")
    ax.set_xlim(0, max(3.0, s.max()))
    ax.legend(loc="upper right", framealpha=0.9, fontsize=8)
    fig.savefig(FIGS / "spacings.pdf")
    plt.close(fig)
    return gaps


def fig_phase():
    """Critical line z P(beta) = 1/4 in the (beta, z) quarter-plane."""
    betas = np.linspace(1.02, 4.0, 200)
    zc = np.array([float(1 / (4 * mp.primezeta(b))) for b in betas])
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(betas, zc, color=BLUE, lw=1.3)
    ax.fill_between(betas, zc, 1e3, color=BLUE, alpha=0.08)
    ax.scatter([SIGMA_STAR], [1.0], color=RED, zorder=5, s=25)
    ax.annotate(r"$(\sigma^{*},1)$", (SIGMA_STAR, 1.0),
                textcoords="offset points", xytext=(8, 6), fontsize=8, color=RED)
    ax.text(3.2, 0.25, "convergent\n(subcritical)", fontsize=8, color=BLUE, ha="center")
    ax.text(1.35, 2.2, "divergent", fontsize=8, color="gray", ha="center")
    ax.set_xlabel(r"inverse temperature $\beta$")
    ax.set_ylabel(r"fugacity $z$")
    ax.set_xlim(1.0, 4.0)
    ax.set_ylim(0, 3.0)
    ax.set_title(r"critical line $z\,P(\beta)=\frac{1}{4}$", fontsize=9)
    fig.savefig(FIGS / "phase.pdf")
    plt.close(fig)


def main():
    # Phase diagram needs no data.
    fig_phase()
    print("wrote phase.pdf")

    census_path = DATA / "zeros_t5000.txt"
    if not census_path.exists() or census_path.stat().st_size == 0:
        print("census file missing/empty; skipping data-driven figures")
        return
    re_parts, im_parts = load_census(census_path)

    # total density D(1) from the curve, else recompute quick
    density_curve = DATA / "density_curve.txt"
    # 4-seed cross-check gives D(1) = 0.14356 +/- 0.00010 (Appendix A.4); the
    # single-run curve file agrees within noise. Quote the cross-check value.
    d1 = 0.1436

    fig_zeros(re_parts, im_parts)
    print("wrote zeros.pdf")
    fig_staircase(im_parts, d1)
    print("wrote staircase.pdf")
    gaps = fig_spacings(im_parts)
    print("wrote spacings.pdf")

    n = len(im_parts)
    print(f"\ncensus summary: {n} zeros, "
          f"Re in ({re_parts.min():.4f}, {re_parts.max():.4f}), "
          f"Im up to {im_parts.max():.1f}")
    print(f"mean spacing {gaps.mean():.4f} (2pi/log2 = {PERIOD:.4f}), "
          f"std {gaps.std():.4f}, CV {gaps.std()/gaps.mean():.4f}")
    print(f"D(1) used for staircase: {d1:.4f}, "
          f"predicted count to Tmax: {d1*im_parts.max():.0f}")


if __name__ == "__main__":
    main()
