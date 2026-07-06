"""Analytic verification for the paper (Appendix A.2-A.3).

1. Critical exponents: P(s)=1/4 (trees), P(s)=1 (words), Otter condition (bushes),
   zeta(s)=2 (Kalmar, for comparison).
2. Numerical check of Z_T(s) = (1 - sqrt(1-4P(s)))/2 against partial sums of g_T(n)/n^s.
3. Numerical check of the Otter equation Z_B(s) = P(s) + (Z_B(s)^2 + Z_B(2s))/2.
4. Complex solutions of P(s) = 1/4 ("tree zeros", branch points of Z_T).
5. Hagedorn-Frautschi-Nahm bootstrap transplanted onto the primon spectrum
   (Boltzmann and Bose counting; see the bootstrap remark and Appendix A).

Requires: pip install mpmath numpy.  Run verify_combinatorics.py first (g_values.txt).
"""
import pathlib
import numpy as np
import mpmath as mp

HERE = pathlib.Path(__file__).parent
mp.mp.dps = 30
P = mp.primezeta

print("=== sanity ===")
print("P(2) =", mp.nstr(P(2), 20), " (known 0.45224742004106549850...)")

print("\n=== critical exponents ===")
s_tree = mp.findroot(lambda s: P(s) - mp.mpf(1) / 4, mp.mpf("2.6"))
s_word = mp.findroot(lambda s: P(s) - 1, mp.mpf("1.4"))
rho_k = mp.findroot(lambda s: mp.zeta(s) - 2, mp.mpf("1.73"))
print("sigma_tree  (P=1/4):", mp.nstr(s_tree, 20))
print("sigma_word  (P=1)  :", mp.nstr(s_word, 20))
print("kalmar rho (zeta=2):", mp.nstr(rho_k, 20))
print("P'(sigma_tree) =", mp.nstr(mp.diff(P, s_tree), 12),
      " cusp amplitude:", mp.nstr(mp.sqrt(-mp.diff(P, s_tree)), 12))
print("P'(sigma_word) =", mp.nstr(mp.diff(P, s_word), 12))

g, u = {}, {}
for line in open(HERE / "g_values.txt"):
    a, b, c = line.split()
    g[int(a)], u[int(a)] = int(b), int(c)

print("\n=== tree zeta identity: partial sum vs (1-sqrt(1-4P))/2 ===")
for s in (3, 4, 6):
    lhs = sum(v * mp.mpf(n) ** (-s) for n, v in g.items())
    rhs = (1 - mp.sqrt(1 - 4 * P(s))) / 2
    print(f"s={s}: partial={mp.nstr(lhs,15)}  closed={mp.nstr(rhs,15)}  diff={mp.nstr(lhs-rhs,3)}")

print("\n=== bush (unordered) zeta via the Otter equation ===")
def Zb(s):
    s = mp.mpf(s)
    ss = [s]
    while ss[-1] < 60:
        ss.append(2 * ss[-1])
    z = P(ss[-1])                      # tail approximation
    for sk in reversed(ss[:-1]):
        z = 1 - mp.sqrt(1 - 2 * P(sk) - z)
    return z

for s in (3, 4):
    lhs = sum(v * mp.mpf(n) ** (-s) for n, v in u.items())
    print(f"s={s}: partial={mp.nstr(lhs,12)}  functional-eq={mp.nstr(Zb(s),12)}  diff={mp.nstr(lhs-Zb(s),3)}")

s_bush = mp.findroot(lambda s: 1 - 2 * P(s) - Zb(2 * s), mp.mpf("1.95"))
print("sigma_bush (1-2P(s)-Z_B(2s)=0):", mp.nstr(s_bush, 20))
print("Z_B at critical point:", mp.nstr(Zb(s_bush + mp.mpf("1e-12")), 10), "(should be ~1)")
print("Hagedorn temperatures:", 1, mp.nstr(1/s_word, 10), mp.nstr(1/s_bush, 10), mp.nstr(1/s_tree, 10))

print("\n=== classical bootstrap on the primon spectrum (Nahm criticality) ===")
TARGET = 2 * mp.log(2) - 1
b_boltz = mp.findroot(lambda s: P(s) - TARGET, mp.mpf("2.15"))
print("beta_H Boltzmann (P = 2log2-1):", mp.nstr(b_boltz, 18))

boot_memo = {}
def G_boot(s):
    """Small branch of the Bose bootstrap G = P + (exp(sum_{k>=1} G(ks)/k) - 1 - G)."""
    s = mp.mpf(s)
    key = mp.nstr(s, 30)
    if key in boot_memo:
        return boot_memo[key]
    if s > 120:
        boot_memo[key] = P(s)
        return boot_memo[key]
    a, Pv = A_boot(s), P(s)
    boot_memo[key] = mp.findroot(lambda g: Pv + mp.e ** (g + a) - 1 - 2 * g, Pv)
    return boot_memo[key]

def A_boot(s):
    s = mp.mpf(s)
    total, k = mp.mpf(0), 2
    while k * s <= 400:
        term = G_boot(k * s) / k
        total += term
        if term < mp.mpf("1e-32"):
            break
        k += 1
    return total

b_bose = mp.findroot(lambda s: P(s) + 2 * A_boot(s) - TARGET, mp.mpf("2.31"))
print("beta_H Bose (P + 2*sum_{k>=2} G(ks)/k = 2log2-1):", mp.nstr(b_bose, 18))
print("G at criticality (1+P)/2:", mp.nstr((1 + P(b_bose)) / 2, 15))

print("\n=== tree zeros: scan of P(s)=1/4 on [1.40,2.70]x[0.3,50] ===")
LIM = 20000
sieve = np.ones(LIM + 1, dtype=bool)
sieve[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i :: i] = False
primes = np.nonzero(sieve)[0].astype(float)
logp = np.log(primes)

t_grid = np.arange(0.3, 50.0, 0.02)
sig_grid = np.arange(1.40, 2.72, 0.02)
phase = np.exp(-1j * np.outer(t_grid, logp))
absF = np.empty((len(sig_grid), len(t_grid)))
for i, sig in enumerate(sig_grid):
    absF[i] = np.abs(phase @ (primes ** (-sig)) - 0.25)

seeds = []
for i in range(1, len(sig_grid) - 1):
    for j in range(1, len(t_grid) - 1):
        v = absF[i, j]
        if v < 0.02 and v <= absF[i - 1 : i + 2, j - 1 : j + 2].min():
            seeds.append(complex(sig_grid[i], t_grid[j]))
print(f"{len(seeds)} candidate seeds")

mp.mp.dps = 20
f = lambda z: P(z) - mp.mpf(1) / 4
roots = []
for sd in seeds:
    try:
        r = mp.findroot(f, mp.mpc(sd), maxsteps=60)
    except Exception:
        continue
    if abs(f(r)) < mp.mpf("1e-15") and 0 < r.imag < 50.5 and 1.0 < r.real < 2.72:
        if not any(abs(r - q) < mp.mpf("1e-6") for q in roots):
            roots.append(r)
roots.sort(key=lambda z: float(z.imag))
period = 2 * mp.pi / mp.log(2)
print("k   sigma                t              t/(2*pi/log 2)")
for k, r in enumerate(roots, 1):
    print(f"{k:2d}  {mp.nstr(r.real,12):20s} {mp.nstr(r.imag,12):14s} {mp.nstr(r.imag/period,6)}")

print("\n=== coarse low-sigma sweep [1.05,1.40) x [0.5,50] ===")
mp.mp.dps = 10
best = (mp.mpf(9), None)
lowseeds = []
for sig100 in range(105, 140, 5):
    sig = mp.mpf(sig100) / 100
    tt = mp.mpf("0.5")
    while tt < 50:
        v = abs(P(mp.mpc(sig, tt)) - mp.mpf(1) / 4)
        if v < best[0]:
            best = (v, mp.mpc(sig, tt))
        if v < mp.mpf("0.05"):
            lowseeds.append(mp.mpc(sig, tt))
        tt += mp.mpf("0.25")
print("min |P-1/4| in band:", mp.nstr(best[0], 5), "at", mp.nstr(best[1], 8))
mp.mp.dps = 20
for sd in lowseeds[:10]:
    try:
        r = mp.findroot(f, sd, maxsteps=60)
        if abs(f(r)) < mp.mpf("1e-15"):
            print("  low-band root:", mp.nstr(r, 12))
    except Exception:
        pass
print("done")
