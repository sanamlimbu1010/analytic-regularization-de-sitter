# Analytic regularization in de Sitter

Companion Mathematica notebooks for Chapter 4 of the MSc thesis *Comparing UV Regularizations for Cosmological Correlators in de Sitter Space* (Sanam Limbu, LMU Munich, 2026). They verify the three-loop pure-necklace difference between analytic and modified analytic regularization: first the k²-derivative (section 4.1), then the closed undifferentiated form (section 4.2).

## Notebooks

| File | Thesis section | Main checks |
|------|----------------|-------------|
| `notebooks/chapter-04-derivative-companion.nb` | 4.1 | Elementary integrals, ZZ term (4.9), kernel Q (4.18), six-Q assembly to boxed derivative (4.21) |
| `notebooks/chapter-04-full-difference-companion.nb` | 4.2 | PZP loop factor (4.23), large-k behaviour, P± and Δ (4.24)–(4.26), derivative vs (4.21), boundary at k → ∞, (4.22) and direct (3.17) numerics |

Each notebook is self-contained. Open in Mathematica 14.3 or later, start a fresh kernel, and evaluate input cells top to bottom. They are saved with outputs so Wolfram Player can display them without re-running.

Exact checks use symbolic simplification under a > 0, b > 0, k > 0 (a = k₁₂, b = k₃₄). Numerical cells state digit agreement reached; that is a regression target, not a certified error bound.

## What is not covered

- Equivalence with other integration orders for [PZP]
- Rigorous bounds behind the large-k outline in the thesis
- The ansatz search that produced the closed form (the result is verified directly)
