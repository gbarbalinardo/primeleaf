"""Tree zeros: census and the Jessen-Tornehave density check.

Implements the numerics behind the paper's sections 5.2-5.3 and Appendix
A.3-A.4: a grid scan of |P(s) - 1/4| for Newton seeds, high-precision
polishing on ``mpmath.primezeta``, and the Monte Carlo estimate of the
mode-dominance probabilities q_p of Theorem 5.3(iii), whose weighted sum
(1/2pi) sum_p log(p) q_p is the predicted zero density.

Requires ``numpy`` and ``mpmath``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np

from primeleaf.primes import primes_up_to

#: Recurrence period of the lightest mode (energy log 2), 2 pi / log 2.
LIGHTEST_MODE_PERIOD = 2 * math.pi / math.log(2)


def _prime_arrays(prime_limit: int) -> Tuple[np.ndarray, np.ndarray]:
    primes = np.array(primes_up_to(prime_limit), dtype=float)
    return primes, np.log(primes)


def scan_seeds(
    sigma_range: Tuple[float, float] = (1.30, 2.72),
    t_range: Tuple[float, float] = (0.3, 500.0),
    *,
    sigma_step: float = 0.02,
    t_step: float = 0.05,
    prime_limit: int = 100_000,
    threshold: float = 0.06,
    block: int = 2000,
) -> List[Tuple[float, float]]:
    """Grid local minima of |P_truncated(s) - 1/4|: Newton seeds for ``polish``."""
    primes, log_primes = _prime_arrays(prime_limit)
    sigma_grid = np.arange(*sigma_range, sigma_step)
    t_grid = np.arange(*t_range, t_step)
    weights = np.array([primes ** (-s) for s in sigma_grid])
    residual = np.empty((len(sigma_grid), len(t_grid)), dtype=np.float32)
    for start in range(0, len(t_grid), block):
        t_block = t_grid[start : start + block]
        phase = np.exp(-1j * np.outer(t_block, log_primes))
        values = phase @ weights.T
        residual[:, start : start + len(t_block)] = np.abs(
            values.T - 0.25
        ).astype(np.float32)
    seeds: List[Tuple[float, float]] = []
    for i in range(1, len(sigma_grid) - 1):
        row = residual[i]
        for j in range(1, len(t_grid) - 1):
            v = row[j]
            if v < threshold and v <= residual[i - 1 : i + 2, j - 1 : j + 2].min():
                seeds.append((float(sigma_grid[i]), float(t_grid[j])))
    return seeds


def polish(
    seeds: Sequence[Tuple[float, float]],
    *,
    dps: int = 18,
    residual: str = "1e-14",
    sigma_box: Tuple[float, float] = (1.0, 2.72),
    t_box: Tuple[float, float] = (0.3, 500.5),
) -> List[mp.mpc]:
    """Newton-polish seeds on the full prime zeta; dedupe, sort by ordinate."""
    zeros: List[mp.mpc] = []
    with mp.workdps(dps):
        tolerance = mp.mpf(residual)

        def f(z):
            return mp.primezeta(z) - mp.mpf(1) / 4

        for sigma, t in seeds:
            try:
                root = mp.findroot(f, mp.mpc(sigma, t), maxsteps=60)
            except Exception:  # findroot: no convergence / division by zero
                continue
            if not abs(f(root)) < tolerance:
                continue
            if not (
                sigma_box[0] < root.real < sigma_box[1]
                and t_box[0] < root.imag < t_box[1]
            ):
                continue
            if any(abs(root - known) < 1e-6 for known in zeros):
                continue
            zeros.append(root)
    zeros.sort(key=lambda z: float(z.imag))
    return zeros


@dataclass
class Census:
    """Result of a tree-zero census (section 5.3, Appendix A.4)."""

    zeros: List[mp.mpc]

    @property
    def real_parts(self) -> List[float]:
        return [float(z.real) for z in self.zeros]

    @property
    def ordinates(self) -> List[float]:
        return [float(z.imag) for z in self.zeros]

    @property
    def spacings(self) -> List[float]:
        t = self.ordinates
        return [b - a for a, b in zip(t, t[1:])]

    def counting(self, height: float) -> int:
        """N(T): number of zeros with ordinate below ``height``."""
        return sum(1 for t in self.ordinates if t < height)

    def save(self, path) -> None:
        lines = (
            f"{mp.nstr(z.real, 15)} {mp.nstr(z.imag, 15)}\n" for z in self.zeros
        )
        Path(path).write_text("".join(lines))

    @staticmethod
    def load(path) -> "Census":
        zeros = []
        for line in Path(path).read_text().splitlines():
            re_part, im_part = line.split()
            zeros.append(mp.mpc(mp.mpf(re_part), mp.mpf(im_part)))
        return Census(zeros)


def census(
    t_max: float = 500.0,
    sigma_range: Tuple[float, float] = (1.30, 2.72),
    *,
    prime_limit: int = 100_000,
    sigma_step: float = 0.02,
    t_step: float = 0.05,
    threshold: float = 0.06,
    dps: int = 18,
) -> Census:
    """Scan and polish all tree zeros in the given box (Appendix A.4)."""
    seeds = scan_seeds(
        sigma_range,
        (0.3, t_max),
        sigma_step=sigma_step,
        t_step=t_step,
        prime_limit=prime_limit,
        threshold=threshold,
    )
    zeros = polish(
        seeds,
        dps=dps,
        sigma_box=(1.0, sigma_range[1]),
        t_box=(0.3, t_max + 0.5),
    )
    return Census(zeros)


def low_band_zeros(
    sigma_range: Tuple[float, float] = (1.05, 1.30),
    t_max: float = 500.0,
    *,
    sigma_step: float = 0.05,
    t_step: float = 0.25,
    scan_dps: int = 8,
    polish_dps: int = 18,
    threshold: float = 0.05,
    max_polish: int = 40,
) -> List[mp.mpc]:
    """Coarse sweep below the main strip (Appendix A.4).  Not exhaustive."""
    seeds: List[Tuple[float, float]] = []
    with mp.workdps(scan_dps):
        quarter = mp.mpf(1) / 4
        sigma = sigma_range[0]
        while sigma < sigma_range[1] - 1e-9:
            t = 0.5
            while t < t_max:
                if abs(mp.primezeta(mp.mpc(sigma, t)) - quarter) < threshold:
                    seeds.append((sigma, t))
                t += t_step
            sigma += sigma_step
    return polish(
        seeds[:max_polish],
        dps=polish_dps,
        sigma_box=(1.0, sigma_range[1] + 0.01),
        t_box=(0.3, t_max + 0.5),
    )


def density_prediction(
    sigma: float = 1.30,
    *,
    prime_limit: int = 20_000,
    samples: int = 400_000,
    track: int = 200,
    block: int = 5_000,
    seed: int = 20260706,
) -> Dict[str, float]:
    """Predicted zero density right of ``sigma``, with an analytic tail.

    Monte Carlo estimate of the dominance probabilities q_p for the first
    ``track`` primes (Theorem 5.3(iii)), a pooled estimate of the local
    density of the random model near 1/4 from the mid-range tracked primes,
    and the analytic tail  sum log(p) pi p^(-2 sigma) f_local  over the
    untracked primes.  Valid down to sigma = 1: the random series converges
    in L^2 for sigma > 1/2 and the density is continuous, so the value at
    sigma = 1 is the total density of tree zeros in Re s > 1.
    """
    rng = np.random.default_rng(seed)
    primes, log_primes = _prime_arrays(prime_limit)
    amplitudes = primes ** (-float(sigma))
    track = min(track, len(primes))
    hits = np.zeros(track)
    done = 0
    while done < samples:
        batch = min(block, samples - done)
        theta = rng.uniform(0.0, 2.0 * np.pi, size=(batch, len(primes)))
        terms = amplitudes * np.exp(1j * theta)
        total = terms.sum(axis=1)
        for k in range(track):
            rest = total - terms[:, k]
            hits[k] += np.count_nonzero(np.abs(rest - 0.25) < amplitudes[k])
        done += batch
    q = hits / samples
    direct = float((log_primes[:track] * q).sum()) / (2 * math.pi)
    window = slice(track // 4, track)
    disk_areas = math.pi * amplitudes[window] ** 2
    f_local = float(q[window].sum() / disk_areas.sum())
    tail_weight = float(
        (log_primes[track:] * math.pi * amplitudes[track:] ** 2).sum()
    )
    tail = f_local * tail_weight / (2 * math.pi)
    return {
        "sigma": float(sigma),
        "direct": direct,
        "local_density": f_local,
        "tail": tail,
        "density": direct + tail,
        "q2": float(q[0]),
        "q3": float(q[1]),
    }


def dominance_probabilities(
    sigma: float = 1.30,
    *,
    prime_limit: int = 10_000,
    samples: int = 400_000,
    track: int = 60,
    block: int = 20_000,
    seed: int = 20260706,
) -> Tuple[Dict[int, float], float]:
    """Monte Carlo q_p of Theorem 5.3(iii) and the predicted zero density.

    q_p is the probability that the p-th mode dominates the coherent rest of
    the random-phase model at height sigma; the predicted density of tree
    zeros to the right of the line Re s = sigma is (1/2pi) sum_p log(p) q_p.
    Returns ``(q, density)``.
    """
    rng = np.random.default_rng(seed)
    primes, _ = _prime_arrays(prime_limit)
    amplitudes = primes ** (-sigma)
    track = min(track, len(primes))
    hits = np.zeros(track)
    done = 0
    while done < samples:
        batch = min(block, samples - done)
        theta = rng.uniform(0.0, 2.0 * np.pi, size=(batch, len(primes)))
        terms = amplitudes * np.exp(1j * theta)
        total = terms.sum(axis=1)
        for k in range(track):
            rest = total - terms[:, k]
            hits[k] += np.count_nonzero(np.abs(rest - 0.25) < amplitudes[k])
        done += batch
    q = {int(primes[k]): float(hits[k] / samples) for k in range(track)}
    density = sum(math.log(p) * v for p, v in q.items()) / (2 * math.pi)
    return q, density
