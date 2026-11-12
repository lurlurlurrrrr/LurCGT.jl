# Part 3b: Numerical irreps constructed from target crystals

Status: implemented. Virtualization supplies the target crystal; numerical
blocks are constructed directly from its strings. The ambient numerical-word
constructor and its dense relation validator have been removed. Sparse
metric-aware generator validation lives in `test/utils/common.jl`.

## Scope and recipes

Construct the B-family vector, G2 seven-dimensional seed, and F4
26-dimensional seed from their target crystals. Keep odd-SO spin seeds
minuscule, including SO(3)'s `(1,)` of dimension 2. Other fundamentals continue
to use recursive product recipes. This direct numerical rule is specific to
these seeds; it is not a general reconstruction algorithm for arbitrary
crystals with weight multiplicities.

Retain the existing recipe:

```julia
struct VirtualRecipe <: FundamentalRecipe
    target::Vector{Int}
    dimension::Int
    ambient::DataType
    operator_words::Vector{Vector{Int}}
end
```

Its ambient data now constructs only the crystal. An entry in
`operator_words[i]` specifies the ordered ambient crystal operations for one
target lowering step. Raising uses the reversed word. Do not request an
ambient numerical `Irep`, apply its numerical matrices, or restrict its metric.
Neither ambient operator-word factorials nor ambient numerical normalization
enter the replacement construction.

Keep `Vector{FundamentalRecipe}`, recursive dependencies, and existing types.
No new recipe fields, catalogue wrappers, or order field are needed.

| Target crystal | Ambient highest label | Operator words |
| --- | --- | --- |
| B_r vector, r≥2 | SO(2r+2), 2F_1 | `[[1,1], …, [r-1,r-1], [r,r+1]]` |
| G2 `(1,0)` | SO(8), `(1,0,1,1)` | `[[1,3,4], [2,2,2]]` |
| F4 `(0,0,0,1)` | E6, `(1,0,0,0,0,1)` | `[[2,2], [4,4], [3,5], [1,6]]` |

Derive the ambient label by adding `target[i]` at every index occurring in
`operator_words[i]`. Retain orbit, repetition, weight-embedding, epsilon/phi,
dimension, and crystal-inverse checks. Respect the existing D3 convention for
SO(5)'s ambient SO(6). Preserve integer-vector crystal words and deterministic
Demazure order, including distinct symbols for equal-weight vertices.

The F4 graph can also be described directly by its 24 short roots and two
zero-weight vertices. This description supplies an independent graph check;
it does not require a new recipe type in this step.

## Numerical construction without weight degeneracy

The B-family vector and G2 seven-dimensional seed have multiplicity-one
weights. Allocate their basis states in the existing target Demazure order,
and group states into `Sz` ranges as before.

For each color i, enumerate maximal target strings

```text
v₀ --i--> v₁ --i--> … --i--> vₙ.
```

Here n counts target crystal edges, not ambient operations or the number of
occurrences of i in a Demazure word. Set

```text
Fᵢ vₖ = (k+1) vₖ₊₁          for 0 ≤ k < n,
⟨vₖ,vₖ⟩ = binomial(n,k)     for 0 ≤ k ≤ n.
```

Start scalar metric blocks at `[1;;]`. A length-two string therefore has
metric entries `1,2,1` and lowering coefficients `1,2`; length-one strings
have metric entries `1,1` and coefficient `1`. Isolated vertices impose no
new norm assignment. Assignments from different nontrivial strings must
agree; report a conflict rather than overwriting an earlier assignment.

These coefficients already encode the target divided-power basis
`vₖ = Fᵢ^k v₀ / k!`. Do not divide the assigned coefficients by another
`div_fac`. Preserve unit highest-state norm, exact integer lowering and Gram
blocks, and the existing primitive inverse-metric representation. Only inverse
metrics may be rational. Raise an error on nonintegral storage or unexpected
multiplicities; do not fall back to ambient numerical construction.

## F4 zero-weight block

For F4 `(0,0,0,1)`, use scalar norm 1 at each of the 24 nonzero short-root
weights. Identify the two zero-weight states by their crystal strings:

```text
x_{α₃} --3--> y₃ --3--> x_{−α₃}
x_{α₄} --4--> y₄ --4--> x_{−α₄}.
```

Keep their actual basis order determined by Demazure traversal. In the semantic
order `(y₃,y₄)`, prescribe the metric

```julia
[2 1; 1 2]
```

Use the same string coefficients as above for the displayed edges, and
coefficient 1 along other short-root crystal edges.

The off-diagonal metric also requires numerical actions that are not crystal
edges. With this positive-overlap convention, the zero-weight blocks
are

```text
F₃ x_{α₃} = y₃,        F₃ (y₃,y₄) = (2,1) x_{−α₃},
F₄ x_{α₄} = y₄,        F₄ (y₃,y₄) = (1,2) x_{−α₄}.
```

The cross coefficients follow from the prescribed metric and requiring the
self-commutator to vanish on the zero-weight space; they must be validated
with the complete matrices. They do not add crystal arrows. Permute these
rows/columns into the actual stored zero-weight order. Do not assume that a
missing crystal edge implies every numerical coefficient is zero in a
repeated-weight block. Keep this exception in the F4 symmetry-specific code.

## Sparse generator validation

Use the metric-aware test utility `lowering_matrix(rep,i,FT)` to obtain a
sparse matrix m_i in an orthonormal basis. Its ordinary adjoint `m_i'` then
represents the metric adjoint. Obtain the ordered stored weight tuples from
`ordered_zvalues(rep)`, and define

```text
Z_a = diag(weight[k][a]),
H_i = m_i' * m_i - m_i * m_i'.
```

For every lowering index i:

1. Check that H_i is diagonal within an explicit numerical tolerance.
2. Solve for coefficients c_i satisfying `H_i = Σ_a c_i[a] Z_a`, and check
   the residual on every basis state. Do not include a constant identity term.
3. Compare c_i across independently constructed minuscule and virtual-crystal
   fundamental seeds of the same symmetry and rank. Use the same node order,
   weight coordinates, normalization, and precision for every comparison.

Fit against the matrix whose rows are ordered z-values. Require enough rank
to determine the coefficients uniquely; a trivial irrep alone cannot establish
a reference. The test utility accepts `reference=coefficients` from a faithful
irrep to validate representations with deficient weight rank, such as the
SO(4) spin fundamentals. Once determined, retain one coefficient vector per
lowering index for that symmetry and test subsequent irreps against it. Product
fundamentals should obey the same coefficients, while independent seeds are
particularly useful for catching inconsistent generator normalizations.

Also compare `Σ_a c_i[a] w[a]` with `getqlabel(S,w)[i]` in the chosen
simple-coroot normalization. This anchors the coefficients when only one
independent seed is available, as for G2 and F4, and prevents consistent but
incorrect rescaling of every seed from passing. Coefficients expressed in
stored z-coordinates need not themselves be integers.

In addition, require these tests for the new seed construction:

1. **Weight shifts.** Every stored nonzero lowering coefficient for index i
   maps weight w to `w - getdz(S,i)`. Check this exactly in the integer blocks;
   weight-block orthonormalization preserves this property.
2. **Mixed commutators.** For every ordered pair `i != j`, check
   `m_i' * m_j - m_j * m_i' ≈ 0`. Include F4's two-dimensional zero-weight
   block and its cross coefficients.
3. **Lowering Serre relations.** Define
   `A[i,j] = getqlabel(S,Tuple(getdz(S,j)))[i]`. Starting with `R = m_j`,
   replace `R` by `m_i * R - R * m_i` exactly `1 - A[i,j]` times and require
   the result to vanish within tolerance. Test every ordered pair `i != j`;
   do not assume the Cartan matrix is symmetric. Nonadjacent nodes require
   one commutator. Raising Serre relations follow by taking adjoints and do
   not need duplicate tests.

Together with the anchored self-commutator check, weight shifts supply the
Cartan-generator relations, while diagonal Cartan operators commute. Thus
these tests cover the defining Lie-algebra relations for the seed matrices.
They supplement dimension, multiplicity, canonical-product, and CG tests;
they do not replace those checks of the intended irrep and basis convention.

Use explicit absolute and relative tolerances scaled by the operand/product
norms for floating-point Cholesky conversion, commutators, and nested Serre
expressions. Compare zero residuals with an explicit norm bound, rather than
using a relative-only comparison against zero. Allow BigFloat checks for
ambiguous residuals. Report the symmetry, label, operator indices, relation,
offending entry or coefficient, and residual on failure.

These checks are required validation tests for the replacement numerical
construction. They replace the previous disabled-relation policy. The old ambient numerical
constructor and its automatic validator are removed.

## Integration and completion

| Family | Direct seeds | Product fundamentals |
| --- | --- | --- |
| SO(3) | Minuscule F_1, dimension 2 | `(2,)`, dimension 3, from two `(1,)` factors. |
| B_r, r≥2 | Crystal vector F_1, dimension 2r+1; minuscule spin F_r, dimension 2^r | F_k from F_(k−1) × F_1, 2≤k≤r−1. |
| G2 | Crystal F_1, dimension 7 | F_2 from F_1 × F_1, dimension 14. |
| F4 | Crystal F_4, dimension 26 | F_1 from F_4 × F_4, dimension 52; F_3 from F_4 × F_4, dimension 273; F_2 from F_3 × F_4, dimension 1274. |

- Replace the ambient numerical path in `virtual.jl`, retaining its crystal
  construction and checks. Put shared target-string numerical helpers in an
  appropriate irrep file and symmetry-specific exceptions in symmetry files.
- Preserve `getNsave_irep` recipe dispatch, deterministic basis order, integer
  storage, recursive products, and unit highest-state norm.
- Test B/G2 binomial norms and coefficients, both F4 zero-weight states and the
  full prescribed Gram block, weight shifts, mixed/self commutators, lowering
  Serre relations, and coefficient agreement
  between B-family vector and minuscule spin seeds.
- Run the same generator validator across all fundamentals of SU(2)–SU(6),
  Sp(2), Sp(4), Sp(6), Sp(8), and SO(6), SO(8), SO(10). For SO(4), establish
  the reference with its vector and check both spin fundamentals. Exceptional
  coverage includes E6 nodes 1,6,2,3,5 (dimensions 27,27,78,351,351) and E7
  nodes 7,1,2,6 (dimensions 56,133,912,1539). Larger exceptional fundamentals
  remain outside routine regression coverage.
- Retain Part 3a validation of incidental components and all product copies
  against canonical construction, comparing every `Irep` field except
  `size_byte`. Retain existing dimension and CG/recoupling tests.
- Use `crystal-seed-v1` for affected B/G2/F4 irreps and dependent CG data,
  in both memory and SQLite. Do not reuse `virtual-divided-v1` numerical data.
  Preserve unchanged minuscule-family caches. No documentation build is needed.

## References

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
chapter 5, supplies the crystal constructions and virtual embeddings.

Benkart, Frenkel, Kang, and Lee,
[Level 1 Perfect Crystals and Path Realizations of Basic Representations at q=0](https://arxiv.org/pdf/math/0507114#page=11),
§3, equations (3.2)–(3.4), gives the root-based little-adjoint crystal, including
the F4 short-root and zero-weight vertices. Translate its twisted-case F4 node
numbering to ours. The numerical normalization and F4 Gram block above are the
implemented normalization policy, checked by the stated validation.
