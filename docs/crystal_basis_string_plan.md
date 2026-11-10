# Part 2: Demazure lowering schedules

Status: implemented following completion of
[part 1](crystal_basis_irrep_plan.md). This document changes no code or caches.

## Goal

Replace the adaptive ascending-color traversal in `get_irep_frommw` with a
fixed sequence of Demazure color stages. Keep `CrystalVertex`, `Irep`, the
vertex/vector dictionary, numerical lowering, and `divcfac`. First discovery
continues to select each vertex's numerical vector and implicit lowering path.
No string-coordinate extraction, new types, or stored graph is required.

## Implementation

1. Hard-code a standard reduced-word formula for each supported symmetry/rank,
   using published conventions and mapping node labels to the local numbering.
   Return colors in chronological lowering order. Validate length and root
   action in tests; no general runtime Weyl-group machinery is needed.
2. Delete the outer queue `while` loop and swap the two `for` loops, replacing
   `1:NZ` with the reduced-word colors. For each color, extend strings from
   **all vertices accumulated before that stage**, including previously visited
   vertices. Use an insertion-ordered vertex vector alongside the existing
   dictionary; freeze its source length at the start of each stage.
3. Walk each source's lowering string. Stop immediately if lowering is undefined
   or the next vertex is already stored, before numerical lowering or `divcfac`.
   Otherwise save the vector and append the vertex in discovery order. The
   `det_getvector` look-ahead is unnecessary.
4. Fill columns within each weight in that explicit discovery order. Retain the
   existing weight ordering and subsequent Gram-matrix/operator construction.

Merely replacing the current inner `1:NZ` loop with a reduced word is not enough:
Demazure stages act on the entire accumulated set, not only the latest queue.
Never mutate vertex words already used as dictionary keys.

The implemented chronological words are:

- A rank r: concatenate `(k, k-1, ..., 1)` for `k=1:r`.
- B/C rank r: repeat `(1, ..., r)` r times.
- D rank r: concatenate `(k, ..., r-2, r-1, r, r-2, ..., k)`
  for `k=1:r-1`, omitting empty ranges.
- G2: `(1, 2, 1, 2, 1, 2)`.

## Ordering and justification

For a written reduced expression ``w_0=s_{i_1}\cdots s_{i_L}``, execute colors
in reverse order: ``i_L,\ldots,i_1``. Starting with ``X_0=\{u_\lambda\}``, each
stage takes the closure under its color:

```math
X_k=\{\widetilde f_{i_{L-k+1}}^a b:
      b\in X_{k-1},\ a\ge0\}\setminus\{0\}.
```

Bump and Schilling, §13.2, equation (13.4) and Theorem 13.5, identify these as
Demazure subsets and guarantee that the final set is the whole normal
highest-weight crystal. In particular, Theorem 13.5(ii) says that the set before
any stage meets each color string in nothing, its highest vertex alone, or the
whole string. Thus a string contributing new vertices starts from its sole
previously present highest vertex; fully present strings contribute nothing.
This justifies retaining first discovery without canonical string extraction.
Fixing the reduced word and source iteration order makes the result reproducible.

Keep the existing numerical lowering, normalization, Gram matrices, and
coefficient solves. Do not add explicit rank checks or a separate independence
gate; let failures in the existing numerical construction propagate. Theorem
13.5(ii) justifies the crystal traversal, rather than independently proving
linear independence of ordinary lowering vectors.

## Validation and handoff

Start with an SU example having a repeated weight. Add small cases to existing
relevant tests for reduced-word length/action, final vertex-set equality with
the old traversal, deterministic first discovery, and generator relations. Include a non-palindromic reduced word to catch
execution-order reversal. Compare different numerical bases through an
intertwining transformation rather than entrywise equality.

Use concrete containers and document new functions. Inspect type stability and
measure allocations/runtime where existing benchmarks permit. Keep independent
per-weight output work separable from the sequential discovery stages.

When activating the implementation, delete the existing representation database
and dependent cached CGT data, then regenerate them under the new basis
convention. This documentation update does not delete databases. Hand the chosen
words and ordering convention to [part 3](crystal_basis_integration_plan.md);
this plan does not require an additional rank-check gate. The cache namespace is
`demazure-v1`, preventing reuse of old files in other installations.

## Reference

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
§13.2, printed p. 174. Supplied PDF:
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

For reduced-word verification, the [SageMath finite Coxeter group reference](https://doc.sagemath.org/html/en/reference/categories/sage/categories/finite_coxeter_groups.html)
provides longest-element and reduced-word operations. Use these as an offline
cross-check when selecting formulas, not as a runtime dependency. Record the
formula source and node-numbering conversion alongside each implementation.
