# Analytic regularization in de Sitter

Companion Mathematica notebooks for the three-loop comparison in Chapter 4 of
the associated MSc thesis.

## Notebooks

- `notebooks/chapter-04-derivative-companion.nb` follows the differentiated
  calculation: the term without a momentum logarithm, the one-logarithm terms,
  the kernel \(\mathcal Q\), the six-\(\mathcal Q\) assembly, and the revised
  boxed derivative.
- `notebooks/chapter-04-full-difference-companion.nb` follows the
  undifferentiated calculation: the PZP prescription, the large-\(k\) boundary
  value, reconstruction through \(P_\pm\), and a direct numerical check.

Each notebook is self-contained. Open it in Mathematica 14.3, start a fresh
kernel, and evaluate the input cells from top to bottom. The last section gives
a compact verification of the main result; preceding sections expose the
intermediate calculations and state the assumptions and evidence status.

The notebooks verify the prescription stated in the thesis. They do not prove
equivalence with other orders of integration for PZP. Numerical precision goals
are regression targets, not certified error bounds.
