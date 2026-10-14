"""Independent numerical checks of Chapter 4 (three-loop difference), in mpmath.

Companion to Chapter 4 of the MSc thesis "Comparing UV Regularizations for
Cosmological Correlators in de Sitter Space" (Sanam Limbu, LMU Munich, 2026).
Equation numbers refer to the thesis; a = k12, b = k34 and k is the exchanged
momentum. Everything below is transcribed from the printed equations of the
thesis; nothing is exported from the Mathematica notebooks. The main checks
(A, B, D) do not use the kernel Q, the assembly, the ansatz search, or the
analytic large-k argument.

  0  Building blocks: the time transforms of 1, log z, log^2 z; the energy
     integrals (loop factors) of (4.23) in closed form against brute-force
     oscillatory quadrature on a grid of z, including k above a and k close to a;
     (4.9) against its s integral;
     the kernel (4.18) against (4.15); the dilogarithm identities used
     (inversion, Landen, the five-term relation (4.28)) and the
     trilogarithm inversion formula.
  A  The derivative (4.21) against the one-dimensional integrals it was
     derived from: the s integral (4.6) and the energy integrals (4.11),
     (4.13), with Euler's constant kept, combined as in (4.2).
  B  The closed form (4.24), (4.25), (4.26) against direct numerical
     integration of the defining integral (3.17): the
     single-energy line with the time integrals first, and [PZP] as in
     (4.23), middle time last. Also against stored values of an earlier,
     separate direct quadrature in Mathematica.
  D  The two properties of the integral used in section 4.2, checked on the
     direct integral itself: its numerical k^2-derivative equals (4.21), and
     it vanishes at large k (k * Delta tends to a constant).
  C  Consistency of the closed form: its k^2-derivative equals (4.21), also
     for k above both momenta; homogeneity of degree -1; symmetry under
     k12 <-> k34; continuity across k12 = k34; stability of the direct
     integration at higher precision.

Checks A, B and C run at fixed points (k below, between and above k12, k34;
k close to k12; k12 close to k34; small and large k) and at seeded random points.
By homogeneity a point is fixed up to scale by r12 = (k12 - k)/(k12 + k) and
r34 = (k34 - k)/(k34 + k) of (4.27), and all positive momenta form the square
-1 < r12, r34 < 1 (physical kinematics is 0 < r12, r34 < 1). The random points
are drawn uniformly from the square |r12|, |r34| < 0.9, with a random overall scale.
Requires mpmath (pip install mpmath). Run: python chapter-04-numerical-checks.py
Runtime: about a minute and a half. A full run reports 132 checks, all passed.
"""
import random

from mpmath import (mp, mpf, mpc, log, pi, polylog, zeta, euler, re, exp, ei, cos,
                    quad, quadosc, inf, diff, sqrt)

mp.dps = 30
Li2 = lambda x: polylog(2, x)
Li3 = lambda x: polylog(3, x)
FAILS = []
COUNT = [0]


def report(label, got, want, tol):
    COUNT[0] += 1
    err = abs(got - want) / max(1, abs(want))
    ok = err < tol
    if not ok:
        FAILS.append(label)
    print(f"  {'ok ' if ok else 'FAIL'} {label}: {mp.nstr(got, 18)} vs {mp.nstr(want, 18)}"
          f" (rel. err {mp.nstr(err, 2)})", flush=True)


def tol(n):
    return mpf(10) ** -n


# ---------------------------------------------------------------- closed forms
def dk_result(a, b, k):
    """d Delta / d k^2, eq. (4.21); a = k12, b = k34."""
    la, lb = log((a + k) / (2 * k)), log((b + k) / (2 * k))
    return pi ** 3 / k * (
        (Li2(-(b - k) / (a + k)) - Li2(-(b - k) / (2 * k)) + la ** 2 / 2) / ((a - k) * (b + k))
        + (Li2(-(a - k) / (b + k)) - Li2(-(a - k) / (2 * k)) + lb ** 2 / 2) / ((a + k) * (b - k))
        + la * lb / ((a - k) * (b - k))
        - (Li2(-(a - k) / (2 * k)) + Li2(-(b - k) / (2 * k)) + pi ** 2 / 6) / ((a + k) * (b + k)))


def delta_closed(a, b, k):
    """Delta L_3, eqs. (4.24), (4.25), (4.26)."""
    la, lb = log((a + k) / (2 * k)), log((b + k) / (2 * k))
    lab = log((a + k) * (b + k) / (4 * k ** 2))
    ra, rb = (a - k) / (a + k), (b - k) / (b + k)
    xa, xb = -(a - k) / (2 * k), -(b - k) / (2 * k)
    Pp = (zeta(3) - Li3(ra) - Li3(rb) - Li3(ra * rb) - 2 * (Li3(xa) + Li3(xb))
          + Li3(-(a - k) / (b + k)) + Li3(-(b - k) / (a + k))
          + lab * (Li2(xa) + Li2(xb)) + lab ** 3 / 6)
    Pm = (2 * (Li3(ra) - Li3(rb)) + 2 * (Li3(xa) - Li3(xb)) - lab * (Li2(xa) - Li2(xb))
          + pi ** 2 / 6 * log((a + k) / (b + k)) - (la ** 3 - lb ** 3) / 3)
    return 2 * pi ** 3 * (Pp / (a + b) + Pm / (a - b))


def zz_closed(a, b, k):
    """[(Z1 - Z4)^2 Pdot], eq. (4.9)."""
    return pi ** 3 / k * (-Li2(-(a - k) / (b + k)) / ((a - k) * (b + k))
                          - Li2(-(b - k) / (a + k)) / ((a + k) * (b - k))
                          + (log((a + k) / (b + k)) ** 2 / 2 + pi ** 2 / 6) / ((a + k) * (b + k)))


def Q_closed(a, b, k):
    """Kernel Q(a, b; k), eq. (4.18)."""
    return pi * ((Li2((b - k) / (a + b)) + log((a + b) * (b + k)) ** 2 / 2 + pi ** 2 / 6) / (b * (a + b))
                 + (Li2((b - a) / (b + k)) - b / (a + b) * (log(a + k) ** 2 - log(b + k) ** 2)) / (b * (a - b)))


# ------------------------------------------------------ time transforms of log^n z
def T(n, a, p):
    """int_0^inf dz e^{-az} log^n(z) cos(pz), n = 0, 1, 2 (Laplace transform at a - ip)."""
    s = mpc(a, -p)
    g = euler + log(s)
    return re([1 / s, -g / s, (g ** 2 + pi ** 2 / 6) / s][n])


def T_direct(n, a, p):
    f = lambda z: exp(-a * z) * log(z) ** n * cos(p * z)
    return quad(f, [0, 1 / p]) + quadosc(f, [1 / p, inf], omega=p)


def line(f, k):
    """Integral over the real p line of an even integrand."""
    return 2 * quad(f, sorted({mpf(0), mpf(1), k, 10 * k, mpf(10), mpf(100)}) + [inf])


# ------------------------------------------------------------------ check A
def derivative_from_integrals(a, b, k):
    """(4.2) from (4.6), (4.11), (4.13) and their images, Euler's constant kept."""
    zz = pi ** 3 / (2 * k) * quad(lambda s: log(s) ** 2 * (
        1 / (a * s + b + k * abs(s - 1)) ** 2 + 1 / ((a + k) * s + b + k) ** 2), [0, 1, 10, inf])

    def C(a, b):  # [Z1 P Pdot], eq. (4.11)
        return pi ** 2 * b * line(
            lambda p: log(p ** 2 + k ** 2) / ((p ** 2 + k ** 2) * (p ** 2 + b ** 2)) * T(1, a, p), k)

    def D(a, b):  # [Pdot Z P], eq. (4.13)
        return -2 * pi ** 2 * b / (k * (a ** 2 - k ** 2)) * line(
            lambda p: log(p ** 2 + k ** 2) / (p ** 2 + b ** 2) * (-a * T(1, k, p) + k * T(1, a, p)), k)

    return -(zz - 2 * (C(a, b) + C(b, a)) + D(a, b) + D(b, a))


# ------------------------------------------------------------------ check B
def F(a, k, z):
    """Loop factor, the energy integral of (4.23) with a = k12 or k34:
    int_R dp a cos(pz) log(p^2+k^2)/(p^2+a^2),
    in closed form. Derived by integrating its k-derivative down from k = infinity;
    checked against brute-force quadrature (F_brute) in section 0."""
    return (pi * exp(-a * z) * log(abs(k ** 2 - a ** 2)) - pi * exp(-a * z) * ei((a - k) * z)
            + pi * exp(a * z) * ei(-(a + k) * z))


def F_brute(a, k, z):
    """The same integral by oscillatory quadrature, after p = u/z so that the
    oscillation has unit period; the peak of width a z at u = 0 is integrated
    separately."""
    f = lambda u: log(u ** 2 / z ** 2 + k ** 2) * cos(u) / (u ** 2 + a ** 2 * z ** 2)
    return 2 * a * z * (quad(f, [0, a * z, 10 * a * z, pi]) + quadosc(f, [pi, inf], omega=1))


def delta_direct(a, b, k):
    """Delta L_3 from (3.17), with the integration orders of chapter 4."""
    P = lambda p: log(p ** 2 + k ** 2)
    zz = pi ** 2 * line(lambda p: P(p) * (T(2, a, p) * T(0, b, p) - 2 * T(1, a, p) * T(1, b, p)
                                          + T(0, a, p) * T(2, b, p)), k)
    zp2 = pi ** 2 * line(lambda p: P(p) ** 2 * (T(1, a, p) * T(0, b, p) + T(0, a, p) * T(1, b, p)), k)
    zs = sorted({mpf(0), 1 / k, mpf(1), mpf(10), mpf(50)}) + [inf]
    pzp = pi * quad(lambda z: 2 * log(z) * F(a, k, z) * F(b, k, z), zs)
    return -(zz - zp2 + pzp)


# ------------------------------------------------------------------ points
FIXED = [("generic", "1.3", "1.7", "0.7"),
         ("unequal", "2.0", "5.0", "0.3"),
         ("k between k12, k34", "0.8", "3.0", "1.5"),
         ("k > k12, k34", "0.8", "1.5", "2.5"),
         ("k12 close to k", "1.000001", "2.2", "1.0"),
         ("k12 close to k34", "1.5", "1.500001", "0.6"),
         ("small k", "3.0", "1.2", "0.05"),
         ("large k", "1.3", "1.7", "100")]
FIXED = [(n, mpf(a), mpf(b), mpf(k)) for n, a, b, k in FIXED]

rng = random.Random(20260930)


def rand_point(i):
    """(r12, r34) uniform in the square |r| < 0.9, k log-uniform in [0.5, 2];
    momenta rounded to 6 significant digits."""
    r12, r34 = rng.uniform(-0.9, 0.9), rng.uniform(-0.9, 0.9)
    k = 0.5 * 4 ** rng.random()
    mom = lambda r: mpf("%.6g" % (k * (1 + r) / (1 - r)))
    return (f"random {i + 1} (r12 = {r12:+.3f}, r34 = {r34:+.3f})", mom(r12), mom(r34), mpf("%.6g" % k))


RANDOM = [rand_point(i) for i in range(12)]
POINTS = FIXED + RANDOM

# ------------------------------------------------------------------ run
print("0  building blocks")
for n in (0, 1, 2):
    for a, p in (("1.3", "0.4"), ("0.7", "5.0")):
        a, p = mpf(a), mpf(p)
        report(f"T_{n}(a={a}, p={p})", T(n, a, p), T_direct(n, a, p), tol(20))
mp.dps = 20
for a, k in (("1.3", "0.7"), ("0.8", "2.5"), ("1.0", "1.000001")):
    for z in ("0.001", "0.03", "0.3", "1", "3", "12", "40"):
        a_, k_, z_ = mpf(a), mpf(k), mpf(z)
        report(f"F(a={a}, k={k}, z={z})", F(a_, k_, z_), F_brute(a_, k_, z_), tol(12))
mp.dps = 30
for n, a, b, k in FIXED[:4]:
    zz_s = pi ** 3 / (2 * k) * quad(lambda s: log(s) ** 2 * (
        1 / (a * s + b + k * abs(s - 1)) ** 2 + 1 / ((a + k) * s + b + k) ** 2), [0, 1, 10, inf])
    report(f"(4.9) {n}", zz_closed(a, b, k), zz_s, tol(20))
for a, b, k in (("1.3", "1.7", "0.7"), ("0.7", "1.7", "0.7"), ("1.3", "0.7", "0.7"), ("0.5", "0.8", "3.0")):
    a, b, k = mpf(a), mpf(b), mpf(k)
    Qq = 2 * quad(lambda p: log(p ** 2 + k ** 2) / (p ** 2 + b ** 2) * re(log(mpc(a, -p)) / mpc(a, -p)),
                  [0, 1, 10, 100, inf])
    report(f"(4.18) a={a}, b={b}, k={k}", Q_closed(a, b, k), Qq, tol(20))
for w in (mpf("-0.3"), mpf("-4.2"), mpf("0.6")):
    report(f"Landen w={w}", Li2(w) + Li2(w / (w - 1)), -log(1 - w) ** 2 / 2, tol(25))
for w in (mpf("-0.3"), mpf("-4.2")):
    report(f"Li2 inversion w={w}", Li2(w) + Li2(1 / w), -pi ** 2 / 6 - log(-w) ** 2 / 2, tol(25))
    report(f"Li3 inversion w={w}", Li3(w), Li3(1 / w) - pi ** 2 / 6 * log(-w) - log(-w) ** 3 / 6, tol(25))
for n, a, b, k in FIXED[:4] + RANDOM[:3]:
    ra, rb = (a - k) / (a + k), (b - k) / (b + k)
    rhs = (Li2(-(a - k) / (b + k)) + Li2(-(b - k) / (a + k)) - Li2(-(a - k) / (2 * k))
           - Li2(-(b - k) / (2 * k)) - log((a + k) / (2 * k)) * log((b + k) / (2 * k)))
    report(f"five-term {n}", Li2(ra * rb), rhs, tol(25))

print("A  derivative (4.21) vs one-dimensional integrals")
for n, a, b, k in POINTS:
    report(n, dk_result(a, b, k), derivative_from_integrals(a, b, k), tol(12))

print("B  closed form (4.24) vs direct integration of (3.17)")
for n, a, b, k in POINTS:
    report(n, delta_closed(a, b, k), delta_direct(a, b, k), tol(12))
print("   vs the earlier Mathematica direct quadrature (18-21 digits)")
stored = [("1.3", "1.7", "0.7", "53.7222810271773134764"),
          ("2", "3", "1", "30.1075194036608193511"),
          ("5", "1.1", "1", "23.8482697266685126321"),
          ("1.3", "1.301", "0.7", "64.8459704100955748964"),
          ("3", "3.5", "0.5", "4.53363179670512925422")]
for a, b, k, v in stored:
    a, b, k = mpf(a), mpf(b), mpf(k)
    report(f"({a}, {b}, {k})", delta_direct(a, b, k), mpf(v), tol(17))

print("D  properties of the integral itself")
for n, a, b, k in FIXED[:5] + RANDOM[:3]:
    report(f"{n}: d/dk^2 of direct integral", diff(lambda q: delta_direct(a, b, sqrt(q)), k ** 2),
           dk_result(a, b, k), tol(12))
a, b = mpf("1.3"), mpf("1.7")
prev = None
for k in (mpf(10), mpf(100), mpf(1000), mpf(10000)):
    d = delta_direct(a, b, k)
    report(f"k = {mp.nstr(k, 1)}: direct vs closed", d, delta_closed(a, b, k), tol(12))
    print(f"       k * Delta = {mp.nstr(k * d, 12)}")
    prev = k * d
report("k * Delta at k = 1e12 (closed) vs 1e4 (direct)", mpf(10) ** 12 * delta_closed(a, b, mpf(10) ** 12),
       prev, tol(3))

print("C  consistency of the closed form")
for n, a, b, k in POINTS:
    report(f"{n}: d/dk^2", diff(lambda q: delta_closed(a, b, sqrt(q)), k ** 2), dk_result(a, b, k), tol(15))
a, b, k = mpf("1.3"), mpf("1.7"), mpf("0.7")
report("homogeneity, factor 3", 3 * delta_closed(3 * a, 3 * b, 3 * k), delta_closed(a, b, k), tol(20))
report("homogeneity, factor 1/7", delta_closed(a / 7, b / 7, k / 7) / 7, delta_closed(a, b, k), tol(20))
report("symmetry k12 <-> k34", delta_closed(b, a, k), delta_closed(a, b, k), tol(20))
report("continuity across k12 = k34", delta_closed(a, a + mpf(10) ** -10, k),
       delta_closed(a, a - mpf(10) ** -10, k), tol(9))
mp.dps = 50
a, b, k = mpf("1.3"), mpf("1.7"), mpf("0.7")
report("direct integration at 50 digits", delta_direct(a, b, k), delta_closed(a, b, k), tol(40))
mp.dps = 30

print(f"\n{COUNT[0]} checks, " + ("ALL PASSED" if not FAILS else f"FAILED: {FAILS}"))
