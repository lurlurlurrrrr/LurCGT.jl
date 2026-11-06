# Part 1: Replace tableau vertices with tensor-word crystal vertices

Status: implementation plan. Preserve the current code structure and replace
`Tableau{S}` as the representation of one crystal vertex. This document changes
no Julia code or cache data.

1. **This part:** one new vertex type and its crystal/tensor operations,
   integrated into the existing traversal for currently supported representations.
2. [Reduced words and string patterns](crystal_basis_string_plan.md): later
   coordinates, basis ordering, and lowering schedules.
3. [Remaining integration](crystal_basis_integration_plan.md): extended seed
   coverage, numerical representations, symmetry contracts, and caches.

## Goal and scope

A vertex of an arbitrary highest-weight crystal should be a point in a tensor
product of fundamental crystals. Each fundamental vertex can in turn be expressed
as a point in an appropriate tensor product of chosen building-block crystals.
Use these embeddings to express a vertex as a tensor word without storing a
Young diagram. This is the general design; the present implementation is limited
to the existing A–D/G2 alphabets and currently supported representation labels.
Spin symbols, spin fundamental words, and numerical spin representations are
deferred to part 3 and are not prerequisites for completing the vertex change.

Keep `get_irep_frommw`'s queue, dictionary of vertices and numerical vectors,
`det_getvector`'s string traversal, and the existing weight-block calculations.
Keep `Crystal_ops{S}` as the operator-data structure. Do not introduce whole-crystal
types, tensor-product container types, graph enumeration as a prerequisite, or a
new numerical representation pipeline.

This part changes vertex storage and the operations on it. It does not introduce
reduced-word schedules or change the chosen numerical basis deliberately. Any
unavoidable change of numerical basis must be identified and handled through
part 3's cache/validation gate rather than treated as a harmless label rename.

## One replacement vertex type

Proposed layout, to be implemented in `src/Symmetry/crystal.jl`:

```julia
struct CrystalVertex{S<:NonabelianSymm}
    word::Vector{Int}
end
```

`word[k]` is an integer encoding one mathematical symbol. For each `S`, the
symbol uniquely determines its weight, elementary raising/lowering transitions,
and epsilon/phi. Reuse the existing symbol encodings and `Crystal_ops` maps;
no separate factor identifiers or factor list are stored. Future alphabets can
use nonoverlapping encodings without changing this type. Mathematical display
symbols need not equal their integer encodings.

The five reasons for this type are:

1. **Why existing structures are insufficient:** `Tableau` makes a shape and
   column-reading word its identity. Tensor factors need not be tableau columns;
   removing the shape requirement makes the point representation independent
   of its tableau construction. Functions alone cannot change those stored semantics.
2. **Invariant:** a vertex is an ordered word of valid local vertices in specified
   building blocks. Operators act on the word and preserve each position's alphabet.
   Reachability from a specified highest word selects the irrep component.
3. **Layout and dispatch:** one concrete integer vector replaces the current
   `shape`/`colread` pair. Only `S` is a type parameter;
   factor count and representation labels stay runtime data. No recursive or
   abstractly typed vertex fields are needed.
4. **Type-stability checks:** inspect lowering, raising, weight lookup, and the
   existing dictionary/queue traversal with `@code_warntype` or JET. Keep
   `Dict{CrystalVertex{S},Vector{RT}}` and `Vector{CrystalVertex{S}}` concrete.
5. **Why this is minimal:** one vertex type replaces one vertex type. The existing
   operator maps, traversal, and numerical storage remain in place; no graph or
   product ownership layer is added.

The user has explicitly requested this vertex-type replacement. Adding other
new types is outside this plan.

Implement equality and hashing by the contents of `word` within symmetry `S`,
not by weight or allocation identity. Distinct vertices at the same weight must
stay distinct. Never mutate the word of a vertex already used as a dictionary
key. Copy it for nonmutating operators; in-place operations are for owned scratch
vertices and must leave them unchanged on failure.

## Fundamental vertices as tensor words

For each currently available fundamental crystal, document its existing embedding

```math
\iota_i:B(\omega_i)\hookrightarrow
C_{i,1}\otimes\cdots\otimes C_{i,\ell_i},
```

where the ``C_{i,k}`` are available building-block crystals of the same type.
Predefine the image of the highest vertex as a word of these symbols. Obtain
the other images by the usual crystal operators; there is no need to pre-enumerate
and store every fundamental vertex before the existing traversal can use it.
The embedding must preserve weights and commute with raising/lowering, including
undefined results. This is a full connected component, not a subset selected
only by weight. If several copies occur, the chosen highest word identifies one.

When all required fundamental crystals are available, for
``\lambda=\sum_i m_i\omega_i`` the highest word is the concatenation of
``m_i`` copies of these fundamental highest words in a fixed order. The component
reached from it is ``B(\lambda)``. A point in the product of fundamental crystals
is represented by concatenating the words of its factor vertices. Associativity
of the tensor rule makes this flattening valid for same-type crystal embeddings;
fundamental block boundaries need not be stored in every vertex. Construction
recipes may retain them as ordinary concrete arrays for validation/debugging.
For λ=0, the word is empty and the weight is zero.

Examples and limits:

- Classical fundamentals already realized inside vector tensor powers can keep
  those realizations. Translate their highest tableaux to words once; tensor
  operations then do not depend on shapes or tableau admissibility predicates.
- In A2, the book convention realizes ``B(\omega_2)`` as the component of
  ``2\otimes1`` in the vector square. This fundamental is a two-letter word;
  it need not be replaced by an opaque graph vertex ID.
- Recursive recipes must terminate at explicit local operator/weight data for
  building blocks. Do not circularly define a missing seed by itself. A known
  seed used as a one-factor crystal is a valid base case.
- A virtual embedding into a different root system is not an ordinary same-type
  tensor embedding. Its grouped/scaled ambient operators must first provide
  validated target-type building-block data. Do not flatten ambient letters into
  an ordinary target-type word and apply the wrong root operators.

Implement the point replacement for A–D and G2 using their existing alphabets.
Preserve the currently supported SO label subset and its vector-derived bootstrap
modules; do not require true spin fundamentals to replace their tableau storage.
Use existing symmetry functions and concrete tuples/arrays; no registry type is
needed. Spin and E/F extensions belong to the later
[seed catalogue](crystal_basis_integration_plan.md#complete-fundamental-construction-catalogue).

## Highest-weight recipes for A–D and G2

For available fundamentals, use a small per-symmetry function taking node k
and returning its highest word. This is a predefined formula or literal word, not a search over
tensor products. The q-label constructor repeats and concatenates these words
where applicable, retaining the existing SO bootstrap recipes described below.
The table uses the book's node numbering and tensor convention; r is the rank.

| Family | Fundamental node k | Highest word | Component dimension |
| --- | --- | --- | --- |
| A_r, `SU{r+1}` | 1≤k≤r | `[k, k-1, ..., 1]` | ``\binom{r+1}{k}`` |
| B_r, `SO{2r+1}` at Lie-algebra level | 1≤k<r | `[k, k-1, ..., 1]` | ``\binom{2r+1}{k}`` |
| C_r, `Sp{2r}` | 1≤k≤r | `[k, k-1, ..., 1]` | ``\binom{2r}{k}-\binom{2r}{k-2}`` |
| D_r, `SO{2r}` at Lie-algebra level | 1≤k≤r−2 | `[k, k-1, ..., 1]` | ``\binom{2r}{k}`` |
| G2 (node 1 short) | k=1 | `[1]` | 7 |
| G2 | k=2 | `[2, 1]` | 14 |

A binomial coefficient with negative lower index is zero. The B/D spin nodes
are deliberately absent from this implementation table.

Retain the existing SO vector-derived highest words as well. In B, the column
`[r, ..., 1]` has weight ``2\omega_r``. In D, `[r-1, ..., 1]` has weight
``\omega_{r-1}+\omega_r``, while `[r, ..., 1]` and `[-r, r-1, ..., 1]`
have weights ``2\omega_r`` and ``2\omega_{r-1}``, respectively. These are
currently useful bootstrap modules, not spin fundamentals. Preserve the existing
`preprocess`/auxiliary-column logic and column order when constructing their words;
only the resulting vertex storage loses its shape field.

Keep existing low-rank handling and accepted labels for B1=A1, C1=A1, D3=A3,
and D2=A1×A1. SO(2) remains outside the simple-family table. Do not add the B1
spin fundamental as part of this refactor; retain its existing vector-derived
construction. No new root-data or label-conversion scheme is required.

### Symbol weights and transformation rules

Retain the current integer vector symbols, including `-j` for barred j.
Only transitions listed below exist; raising maps are their inverses. Derive
local epsilon/phi from the complete monochromatic strings, as in the existing
`maximal_apply_num`. In B and G2, symbol `0` is a real vertex; use missing map
entries / `nothing` for undefined action, never the integer zero as a sentinel.

| Family | Local lowering transitions |
| --- | --- |
| A_r | ``f_i(i)=i+1`` for 1≤i≤r on alphabet 1,…,r+1. |
| B_r vector | For i<r, ``i\mapsto i+1`` and ``\overline{i+1}\mapsto\bar i``; at r, ``r\mapsto0\mapsto\bar r``. |
| C_r vector | For i<r, the same paired transitions; at r, ``r\mapsto\bar r``. |
| D_r vector | For i<r, the same paired transitions; at r, ``r-1\mapsto\bar r`` and ``r\mapsto\overline{r-1}``. |
| G2 | ``f_1:1\mapsto2,\ 3\mapsto0\mapsto\bar3,\ \bar2\mapsto\bar1``; ``f_2:2\mapsto3,\ \bar3\mapsto\bar2``. |

Reuse the existing vector-symbol weights and maps in `SU.jl`, `SO.jl`, `Sp.jl`,
and `G2.jl`. In particular, the G2 highest words agree with the existing
`mw_column` values; its second word has z-weight `(0,2)`, matching q-label `(0,1)`.

Validate highest-word weights against the existing label conversions. Preserve
SO tableau label preprocessing for its supported vector-derived modules; changing
spin-label acceptance or repairing unrelated conversion logic is outside this
point-type refactor.

For a currently supported q-label, concatenate the prescribed words in a documented fixed
order. First match existing column order wherever possible so this vertex change
does not gratuitously change numerical traversal. Test every highest-word recipe
for weight, annihilation by all raising operators, and expected component size
at small ranks. Do not require materializing every rank's fundamental crystals
merely to make the recipes available.

## Operators and compatibility with current code

Extend the current operations to `CrystalVertex{S}`:

| Existing operation | Minimal adaptation |
| --- | --- |
| `mw_tableau_qlabel` | Add a highest-vertex constructor using existing fundamental and SO bootstrap words. Retain label preprocessing and column order. |
| `get_crystal_ops`, `get_eops_std`, `maximal_apply_num` | Reuse the current alphabet maps. Validate inverse transitions and string lengths; reject collisions or inconsistent reverse edges. |
| `apply_lowering!`, `apply_lowering` | Apply the tensor rule to `word`; preserve the existing `Bool` / vertex-or-`nothing` result conventions and failure behavior. |
| Raising and epsilon/phi access | Add the corresponding vertex operations using the same alphabet maps and tensor convention. |
| `get_weight` | Sum the existing symbol weights in the current z-coordinate convention; keep the current alphabet-to-weight lookup. |
| `det_getvector` | Change its vertex argument and dictionary key type; retain its current search along a lowering string. |
| `get_irep_frommw` | Replace the highest-tableau construction, queue element type, and dictionary key type. Retain traversal order, numerical vector propagation, normalization, and weight-block assembly. |

Do not impose a new internal weight-coordinate system in this part. Check
q-label/highest-weight agreement through existing adapters. Broader root-data
and public spin-label changes remain in part 3.

Adopt Bump–Schilling §2.3's left-to-right tensor convention explicitly:

```math
f_i(x\otimes y)=
\begin{cases}
f_i(x)\otimes y,&\phi_i(y)\leq\epsilon_i(x),\\
x\otimes f_i(y),&\phi_i(y)>\epsilon_i(x),
\end{cases}
\qquad
e_i(x\otimes y)=
\begin{cases}
e_i(x)\otimes y,&\phi_i(y)<\epsilon_i(x),\\
x\otimes e_i(y),&\phi_i(y)\geq\epsilon_i(x).
\end{cases}
```

Weights add, and

```math
\phi_i(x\otimes y)=\phi_i(x)+\max(0,\phi_i(y)-\epsilon_i(x)),\qquad
\epsilon_i(x\otimes y)=\epsilon_i(y)+\max(0,\epsilon_i(x)-\phi_i(y)).
```

For multiple factors use the associative rule, or the equivalent §2.4 signature:
each factor contributes ``-^{\phi_i}+^{\epsilon_i}``; cancel adjacent `+-`
pairs, apply f to the rightmost remaining `-`, and e to the leftmost remaining
`+`. Use integer counts and reusable workspace rather than allocating a string
of signs for each call. Change exactly one factor vertex on success.

Do not copy the existing column-reading order implicitly: the book uses the
opposite convention from Kashiwara's original papers. Establish an explicit
factor-order adapter when comparing with existing tableaux. For A1 with
``f(u)=v``, the book convention gives the chain
``u\otimes u\to u\otimes v\to v\otimes v`` and isolated ``v\otimes u``.
This small example fixes both tie cases and factor order.

The existing implementation selects an index in `colread`; check its convention
against the formulas above before translating it to `word`. If it stores factors
in the reverse order, adapt the representation explicitly while preserving the
current traversal's meaning. Remove tableau admissibility checks from the new
vertex operations only where the proven component embedding supplies closure.
Local inverse-edge and string-length checks alone do not prove a seed is normal.

## Implementation and validation sequence

1. Add `CrystalVertex{S}`, content equality/hash, and the local alphabet/weight
   adapters. Add concise docstrings for every new named function.
2. Translate the existing A–D/G2 highest-vertex construction to words, including
   SO bootstrap/auxiliary columns. At small ranks compare reachable vertices,
   weights, and transitions with the existing tableaux; include the 14-vertex
   G2 component. No spin alphabet or spin test fixture is needed.
3. Adapt lowering and add raising/epsilon/phi. Test the A1 convention example,
   both tie cases, inverse operators, empty/single-factor words, and failure
   without mutation. Check grouped fundamental operations agree with flattened
   words; include A2 ``B(\omega_1+\omega_2)`` with its repeated zero weight.
4. Replace the vertex type in the existing `det_getvector` and `get_irep_frommw`
   paths. Do not split graph discovery from vector construction in this step.
   Compare dimensions, exact block ranks, metrics and numerical generator action;
   compare through a basis map if the old and new column order differ.
5. Verify supported SO labels retain their existing behavior and unsupported
   spin requests remain unsupported. Keep spin symbol maps, fundamental recipes,
   numerical seeds, and public-label extensions entirely in the later work.

Treat the behavior-preserving vertex replacement as an internal refactor: do not
add tests or fixtures by default, and run only the prescribed DMRG test runner.
The checks above specify invariants to preserve during review. If implementation
requires new behavior, separate it from this refactor and follow the normal
relevant-test policy. Inspect type stability and measure allocation/runtime
changes against the current tableau traversal.
Preserve opportunities for independent per-weight numerical work without
rewriting the existing traversal solely to parallelize it. Apply the
[shared implementation policy](crystal_basis_integration_plan.md#shared-implementation-policy).

Completion means the existing construction can use word-only vertices in place
of tableaux for currently supported A–D/G2 representations, existing highest
words and tensor embeddings are preserved, and consumers use the
same weight/raising/lowering operations. Part 2 consumes these vertex operations
directly; it does not require owned whole-crystal graphs or precomputed edge tables.

## References and provenance

Primary reference: Daniel Bump and Anne Schilling, *Crystal Bases:
Representations and Combinatorics* (2017), the supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

- §2.2, Definition 2.13 and equation (2.6): crystal axioms and string lengths.
- §§2.3–2.4, printed pp. 18–23: tensor convention, associativity, and signatures.
- §5.4, Proposition 5.11: classical fundamental words and the distinction
  between spin fundamentals and vector-derived orthogonal bootstrap modules.
- §5.7, Theorem 5.20 and its proof, printed pp. 72–73: normal crystals and the
  highest component in tensor products of fundamental crystals.
- Chapters 3 and 6: tableau realizations retained as seed builders/references.

These passages were consulted in the local PDF's extracted text for this rewrite.
The vertex type and storage choices are engineering proposals, not prescriptions from
the book. This documentation-only change ran no numerical tests or doc builds.
