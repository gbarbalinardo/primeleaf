"""Tree-zero census to height 500 and Jessen-Tornehave density check (paper Sec. 5.3).

1. Grid scan of |P(s) - 1/4| (primes <= 10^5) over [1.30, 2.72] x [0.3, 500],
   Newton-polished on mpmath.primezeta; writes zeros_t500.txt.
2. Coarse low-sigma sweep [1.05, 1.30) x [0.5, 500] for stray zeros.
3. Monte Carlo estimate of the Jessen-Tornehave density at sigma_1 = 1.30 via the
   mode-dominance formula  D = (1/2pi) * sum_p log(p) * q_p(sigma_1),
   q_p = Prob(|X - a_p e^{i theta_p} - 1/4| < a_p),  X = sum_p a_p e^{i theta_p}.

Requires: pip install mpmath numpy.
"""
import pathlib, math
import numpy as np
import mpmath as mp

HERE = pathlib.Path(__file__).parent

# ---------- primes ----------
LIM = 100000
sieve = np.ones(LIM + 1, dtype=bool)
sieve[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i :: i] = False
primes = np.nonzero(sieve)[0].astype(float)
logp = np.log(primes)
print(f"{len(primes)} primes up to {LIM}", flush=True)

# ---------- grid scan ----------
SIG0, SIG1, DS = 1.30, 2.72, 0.02
T0, T1, DT = 0.3, 500.0, 0.05
sig_grid = np.arange(SIG0, SIG1, DS)
t_all = np.arange(T0, T1, DT)
W = np.array([primes ** (-s) for s in sig_grid])          # (S, P)
absF = np.empty((len(sig_grid), len(t_all)), dtype=np.float32)
B = 2000
for b0 in range(0, len(t_all), B):
    tb = t_all[b0 : b0 + B]
    phase = np.exp(-1j * np.outer(tb, logp))              # (B, P)
    vals = phase @ W.T                                     # (B, S)
    absF[:, b0 : b0 + len(tb)] = np.abs(vals.T - 0.25).astype(np.float32)
    print(f"grid block t={tb[0]:.1f}..{tb[-1]:.1f}", flush=True)

seeds = []
for i in range(1, len(sig_grid) - 1):
    row = absF[i]
    for j in range(1, len(t_all) - 1):
        v = row[j]
        if v < 0.06 and v <= absF[i - 1 : i + 2, j - 1 : j + 2].min():
            seeds.append((sig_grid[i], t_all[j]))
print(f"{len(seeds)} candidate seeds", flush=True)

mp.mp.dps = 18
f = lambda z: mp.primezeta(z) - mp.mpf(1) / 4
roots = []
for sg, tt in seeds:
    try:
        r = mp.findroot(f, mp.mpc(sg, tt), maxsteps=60)
    except Exception:
        continue
    if abs(f(r)) < mp.mpf("1e-14") and 0.3 < r.imag < 500.5 and 1.0 < r.real < 2.72:
        if not any(abs(r - q) < 1e-6 for q in roots):
            roots.append(r)
roots.sort(key=lambda z: float(z.imag))
with open(HERE / "zeros_t500.txt", "w") as fh:
    for r in roots:
        fh.write(f"{mp.nstr(r.real, 15)} {mp.nstr(r.imag, 15)}\n")

n = len(roots)
res = [float(r.real) for r in roots]
ims = [float(r.imag) for r in roots]
period = 2 * math.pi / math.log(2)
print(f"\n{n} zeros with 0 < t < 500, Re s in ({min(res):.4f}, {max(res):.4f})", flush=True)
print("first five:", [(round(a, 5), round(b, 5)) for a, b in zip(res[:5], ims[:5])])
sp = np.diff(ims)
print(f"mean spacing {np.mean(sp):.4f}  (2pi/log2 = {period:.4f}),  std {np.std(sp):.4f},  "
      f"min {min(sp):.3f}, max {max(sp):.3f}")
for Tc in (100, 200, 300, 400, 500):
    c = sum(1 for t in ims if t < Tc)
    print(f"N({Tc}) = {c:3d}    T*log2/2pi = {Tc/period:.2f}")

# ---------- low-sigma band ----------
print("\nlow band [1.05,1.30) x [0.5,500], coarse mpmath sweep", flush=True)
mp.mp.dps = 8
lowseeds = []
best = 9.0
sig = 1.05
while sig < 1.295:
    tt = 0.5
    while tt < 500:
        v = abs(mp.primezeta(mp.mpc(sig, tt)) - 0.25)
        if v < best:
            best = float(v)
        if v < 0.05:
            lowseeds.append((sig, tt))
        tt += 0.25
    sig += 0.05
print(f"min |P-1/4| in band: {best:.4f};  {len(lowseeds)} sub-0.05 dips", flush=True)
mp.mp.dps = 18
lowfound = []
for sg, tt in lowseeds[:40]:
    try:
        r = mp.findroot(f, mp.mpc(sg, tt), maxsteps=60)
        if abs(f(r)) < mp.mpf("1e-14") and r.real < 1.31 and 0.3 < r.imag < 500.5:
            if not any(abs(r - q) < 1e-6 for q in lowfound):
                lowfound.append(r)
                print("  LOW-SIGMA ZERO:", mp.nstr(r, 12))
    except Exception:
        pass
if not lowfound:
    print("  none polished into the band Re s < 1.30")

# ---------- Monte Carlo Jessen-Tornehave density at sigma = 1.30 ----------
print("\nMonte Carlo dominance probabilities at sigma = 1.30", flush=True)
rng = np.random.default_rng(20260706)
pm = primes[primes <= 10000]                      # 1229 primes
amp = pm ** (-1.30)
NS, BLK = 400000, 20000
KTRACK = 60                                        # track q_p for first 60 primes
hits = np.zeros(KTRACK)
for _ in range(NS // BLK):
    theta = rng.uniform(0, 2 * np.pi, size=(BLK, len(pm)))
    terms = amp * np.exp(1j * theta)
    X = terms.sum(axis=1)
    for k in range(KTRACK):
        rest = X - terms[:, k]
        hits[k] += np.count_nonzero(np.abs(rest - 0.25) < amp[k])
q = hits / NS
D = sum(math.log(pm[k]) * q[k] for k in range(KTRACK)) / (2 * math.pi)
print("q_p for p=2,3,5,7,11,13:", [round(float(x), 4) for x in q[:6]])
print(f"predicted density D(1.30) = {D:.5f} per unit t  ->  {500*D:.1f} zeros to t=500")
print(f"observed: {n};  one-per-period would be {500/period:.1f}")
print("done", flush=True)
