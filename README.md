# Analytic regularization in de Sitter: Chapter 4 companion

Mathematica notebooks and a Python script that accompany Chapter 4 of the MSc thesis
*Comparing UV Regularizations for Cosmological Correlators in de Sitter Space*
(Sanam Limbu, LMU Munich, 2026). They check, by computer algebra and by independent
numerical integration, the closed-form results of that chapter. You do not need the
thesis to run them, but you need it to read what is being checked.

## Context

Cosmological correlators are the basic observables of cosmology. The thesis studies
them for a conformally coupled scalar with a $\phi^4$ interaction in four-dimensional
de Sitter space. Loop integrals over the conformal times of the vertices diverge in the
ultraviolet, so a UV regularization is needed. Two de Sitter-invariant choices are compared:

- **Analytic regularization** (AR) deforms the conformal dimensions in the shadow action.
- **Modified analytic regularization** (MAR) moves the regulator's dependence on the
  conformal times of the internal vertices to those of the external vertices. This makes it
  possible to resum necklace diagrams to all orders.

A *necklace* diagram of the four-point function is a chain of one-loop bubbles between the
vertex that carries $k_1,k_2$ and the vertex that carries $k_3,k_4$. For necklaces the two schemes
agree at one and two loops and first differ at three loops. Chapter 4 computes that
three-loop difference $\Delta L_3$ in closed form, by first computing its derivative with
respect to $k^2$ and then reconstructing $\Delta L_3$ itself.

**Notation.** $k=|\vec k_1+\vec k_2|=|\vec k_3+\vec k_4|$ is the exchanged momentum,
$k_{12}=|\vec k_1|+|\vec k_2|$ and $k_{34}=|\vec k_3|+|\vec k_4|$. Kinematics is
$k>0$, $k_{12}>k$, $k_{34}>k$; the closed forms hold for all positive $k_{12},k_{34},k$.
In the notebooks and the script, `a` $=k_{12}$, `b` $=k_{34}$ and `k` $=k$. Equation numbers
are those of the submitted thesis.

### The results that are checked

With $r_{12}=\frac{k_{12}-k}{k_{12}+k}$ and $r_{34}=\frac{k_{34}-k}{k_{34}+k}$, the difference is

```math
\Delta L_3=2\pi^3\Big[\frac{P_+}{k_{12}+k_{34}}+\frac{P_-}{k_{12}-k_{34}}\Big]\qquad(4.24)
```

with trilogarithms in $P_\pm$ (eqs. (4.25) and (4.26) of the thesis). Its derivative is
(eq. (4.21))

```math
\begin{aligned}
\partial_{k^2}\Delta L_3=\frac{\pi^3}{k}\Big[
&\frac{\mathrm{Li}_2\big(-\frac{k_{34}-k}{k_{12}+k}\big)-\mathrm{Li}_2\big(-\frac{k_{34}-k}{2k}\big)+\frac12\log^2\frac{k_{12}+k}{2k}}{(k_{12}-k)(k_{34}+k)}\\
&+\frac{\mathrm{Li}_2\big(-\frac{k_{12}-k}{k_{34}+k}\big)-\mathrm{Li}_2\big(-\frac{k_{12}-k}{2k}\big)+\frac12\log^2\frac{k_{34}+k}{2k}}{(k_{12}+k)(k_{34}-k)}\\
&+\frac{\log\frac{k_{12}+k}{2k}\,\log\frac{k_{34}+k}{2k}}{(k_{12}-k)(k_{34}-k)}\\
&-\frac{\mathrm{Li}_2\big(-\frac{k_{12}-k}{2k}\big)+\mathrm{Li}_2\big(-\frac{k_{34}-k}{2k}\big)+\frac{\pi^2}{6}}{(k_{12}+k)(k_{34}+k)}\Big]
\end{aligned}
```
The step from the derivative to $\Delta L_3$ uses two properties of the integral: its
$k^2$-derivative is (4.21), and it vanishes as $k\to\infty$ at fixed $k_{12},k_{34}$.
$\Delta L_3$ is the integral (3.17) of the thesis, with the term $[PZP]$ taken in the
integration order of eq. (4.23).

## Files

| File | Thesis section | What it checks |
|------|----------------|----------------|
| `chapter-04-derivative-companion.nb` | 4.1 | $[(Z_1-Z_4)^2\dot P]$ in closed form (4.9) against its $s$ integral (4.6); the kernel $\mathcal Q$: primitive (4.17), closed form (4.18), derivative (4.16), and (4.18) against quadrature of (4.15); the exact assembly of the derivative from six values of $\mathcal Q$ and (4.9) through (4.14), with the ten dilogarithm arguments, the Landen pairs (4.19)-(4.20) and the four remaining dilogarithms, reproducing (4.21); (4.21) against the one-dimensional integrals (4.6), (4.11), (4.13) |
| `chapter-04-full-difference-companion.nb` | 4.2 | the closed form (4.24)-(4.26) typed as printed; Landen (4.19) for $r_{12},r_{34}$ and the five-term relation (4.28); the exact $k$-derivative of (4.24), which equals $2k$ times (4.21); the limits of the trilogarithm arguments ($-1$, $\frac12$, $1$, $1$) and $P_\pm\to0$; (4.24) against the integral (4.22); the loop factor of (4.23) against oscillatory quadrature; direct integration of (3.17) with $[PZP]$ as in (4.23) against (4.24) |
| `chapter-04-numerical-checks.py` | 4.1, 4.2 | independent numerical checks (132) with mpmath at 30 digits: building blocks, (4.21) against the integrals (4.6), (4.11), (4.13); (4.24)-(4.26) against direct integration of (3.17) with $[PZP]$ as in (4.23), at fixed and random points of $-1<r_{12},r_{34}<1$; the $k^2$-derivative of the direct integral against (4.21) and its fall-off like $1/k$; consistency of the closed form |
| `requirements.txt` | | `mpmath>=1.3` |

The script is transcribed from the printed equations of the thesis. It does not use the
notebooks, the kernel $\mathcal Q$, the assembly, or the ansatz search for its main checks.

## How to run

**Notebooks.** Use Mathematica 14.3 or later. Open a notebook, start a fresh kernel and
evaluate all input cells from top to bottom. The first input cell clears the `Global` context. Each
notebook takes under a minute. They are saved with their outputs, so the free
[Wolfram Player](https://www.wolfram.com/player/) shows the results without Mathematica.
Each section ends with a `PASS` line and the last cell is a ledger of all checks, with the kind
of check (exact or numerical) and the residual or the number of digits of agreement. A failed
check aborts the evaluation.

**Script.** Python 3 and mpmath:

```sh
pip install -r requirements.txt
python chapter-04-numerical-checks.py
```

It takes about 2 minutes and ends with `132 checks, ALL PASSED`.

## What is exact and what is numerical

- **Exact** (notebooks): identities that reduce to $0$ symbolically for generic positive
  $k_{12},k_{34},k$, namely the $k$-derivative of (4.18), the six-$\mathcal Q$ assembly to (4.21), the
  Landen pairs, the five-term relation (4.28) (both sides vanish at $k_{34}=k$ and the
  derivatives agree), the derivative of (4.24) against $2k$ times (4.21), and the special values at large $k$.
  Logarithms with positive arguments are split with `PowerExpand`, so these checks assume the
  positivity of the arguments, which holds for positive momenta.
- **Numerical** (notebooks and script): agreement between a closed form and a numerical
  integral or an identity evaluated at points, including points with $k>k_{12},k_{34}$. The
  notebooks record the digits of agreement reached (12 to 40 digits), the script a relative error
  below its tolerance at 30 digits working precision. A precision goal is a target for the
  quadrature, not a certified error bound.

## What is not covered

- Equivalence with other integration orders for $[PZP]$. The results hold for the order of
  eq. (4.23); some other orders diverge.
- The bounds behind the large-$k$ limit of $\Delta L_3$. The thesis argues the limit in outline.
  The notebook checks only the limits of the closed form, and the script checks numerically that the direct
  integral falls off like $1/k$.
- The ansatz search that found the closed form (4.24). It is not part of the argument, since
  the result is checked directly against (3.17).
