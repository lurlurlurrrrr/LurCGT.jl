# Part 1: Replace tableau vertices with tensor-word crystal vertices

Status: implementation plan. Preserve the current code structure and replace
`Tableau{S}` as the representation of one crystal vertex. This document changes
no Julia code or cache data.

1. **This part:** one new vertex type, fundamental-vertex tensor embeddings,
   and the usual crystal/tensor operators within the existing traversal.
2. [Reduced words and string patterns](crystal_basis_string_plan.md): later
   coordinates, basis ordering, and lowering schedules.
3. [Remaining integration](crystal_basis_integration_plan.md): extended seed
   coverage, numerical representations, symmetry contracts, and caches.

## Goal and scope

A vertex of an arbitrary highest-weight crystal should be a point in a tensor
product of fundamental crystals. Each fundamental vertex can in turn be expressed
as a point in an appropriate tensor product of chosen building-block crystals.
Use these embeddings to express a vertex as a tensor word, without requiring a
Young diagram or spin half-column in its storage.

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
    factors::Vector{Int}
    word::Vector{Int}
end
```

`factors[k]` identifies the building-block crystal used at tensor position k;
`word[k]` identifies its current vertex. Both vectors have the same length.
Use stable, disjoint integer alphabets for the building blocks of a given `S`,
so the existing `Crystal_ops` integer-keyed transition and string-length maps
can address their union without changing their field types. Factor identifiers
also identify the ambient product for equality and validation. A fundamental
label is not necessarily a building-block identifier.

The five reasons for this type are:

1. **Why existing structures are insufficient:** `Tableau` makes a shape and
   column-reading word its identity. Tensor factors need not be tableau columns;
   assigning artificial shapes to spin or exceptional vertices preserves the
   unwanted restriction. Functions alone cannot change those stored semantics.
2. **Invariant:** a vertex is an ordered word of valid local vertices in specified
   building blocks. Operators preserve the factor list and act on the word.
   Reachability from a specified highest word selects the irrep component.
3. **Layout and dispatch:** the two concrete integer vectors mirror the size and
   simplicity of the current `shape`/`colread` layout. Only `S` is a type parameter;
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

Implement equality and hashing by contents of both `factors` and `word`, not by
weight or allocation identity. Distinct vertices at the same weight must stay
distinct. Never mutate either vector of a vertex already used as a dictionary
key. Copy the word for nonmutating operators; a factor list may be shared only
under a documented read-only ownership rule. In-place operations are for owned
scratch vertices and must leave them unchanged on failure.

## Fundamental vertices as tensor words

For each required fundamental crystal choose and document an embedding

```math
\iota_i:B(\omega_i)\hookrightarrow
C_{i,1}\otimes\cdots\otimes C_{i,\ell_i},
```

where the ``C_{i,k}`` are available building-block crystals of the same type.
Store the ordered factor identifiers and the image of the highest vertex. Obtain
the other images by the usual crystal operators; there is no need to pre-enumerate
and store every fundamental vertex before the existing traversal can use it.
The embedding must preserve weights and commute with raising/lowering, including
undefined results. This is a full connected component, not a subset selected
only by weight. If several copies occur, the chosen highest word identifies one.

For ``\lambda=\sum_i m_i\omega_i``, the highest word is the concatenation of
``m_i`` copies of these fundamental highest words in a fixed order. The component
reached from it is ``B(\lambda)``. A point in the product of fundamental crystals
is represented by concatenating the words of its factor vertices. Associativity
of the tensor rule makes this flattening valid for same-type crystal embeddings;
fundamental block boundaries need not be stored in every vertex. Construction
recipes may retain them as ordinary concrete arrays for validation/debugging.
For λ=0, both vectors are empty and the weight is zero.

Examples and limits:

- Classical fundamentals already realized inside vector tensor powers can keep
  those realizations. Translate their highest tableaux to words once; tensor
  operations then do not depend on shapes or tableau admissibility predicates.
- In A2, the book convention realizes ``B(\omega_2)`` as the component of
  ``2\otimes1`` in the vector square. This fundamental is a two-letter word;
  it need not be replaced by an opaque graph vertex ID.
- Spin crystals require spin-capable building blocks: vector tensor powers cannot
  produce their weight-lattice coset. A directly known spin seed is a local
  alphabet, with each seed vertex a one-factor word. Higher products use exactly
  the same vertex type and tensor rule; no half-column field is needed.
- Recursive recipes must terminate at explicit local operator/weight data for
  building blocks. Do not circularly define a missing seed by itself. A known
  seed used as a one-factor crystal is a valid base case.
- A virtual embedding into a different root system is not an ordinary same-type
  tensor embedding. Its grouped/scaled ambient operators must first provide
  validated target-type building-block data. Do not flatten ambient letters into
  an ordinary target-type word and apply the wrong root operators.

Use existing symmetry functions and concrete tuples/arrays for these recipes;
no registry type is needed. The extended
[seed catalogue](crystal_basis_integration_plan.md#complete-fundamental-construction-catalogue)
is separate from this small vertex/traversal change. Missing recipes should
report the unsupported fundamental node, not silently substitute another module.

## Operators and compatibility with current code

Extend the current operations to `CrystalVertex{S}`:

| Existing operation | Minimal adaptation |
| --- | --- |
| `mw_tableau_qlabel` | Add a highest-vertex constructor that assembles the embedded fundamental highest words. Retain existing label conversion conventions. |
| `get_crystal_ops`, `get_eops_std`, `maximal_apply_num` | Reuse the current maps for the disjoint building-block alphabets. Validate inverse transitions and string lengths; reject collisions or inconsistent reverse edges. |
| `apply_lowering!`, `apply_lowering` | Apply the tensor rule to `word`; preserve the existing `Bool` / vertex-or-`nothing` result conventions and failure behavior. |
| Raising and epsilon/phi access | Add the corresponding vertex operations using the same alphabet maps and tensor convention. |
| `get_weight` | Sum local building-block weights in the existing z-coordinate convention. The current defining-vector-only lookup must accept the selected alphabets. |
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
2. Translate existing highest-tableau constructors into fundamental highest-word
   recipes for small supported cases. Compare the complete reachable vertices,
   weights, and colored transitions under the explicit embedding.
3. Adapt lowering and add raising/epsilon/phi. Test the A1 convention example,
   both tie cases, inverse operators, empty/single-factor words, and failure
   without mutation. Check grouped fundamental operations agree with flattened
   words; include A2 ``B(\omega_1+\omega_2)`` with its repeated zero weight.
4. Replace the vertex type in the existing `det_getvector` and `get_irep_frommw`
   paths. Do not split graph discovery from vector construction in this step.
   Compare dimensions, exact block ranks, metrics and numerical generator action;
   compare through a basis map if the old and new column order differ.
5. Exercise a small spin-capable alphabet and existing G2 data through the same
   vertex operations. New public spin support and exceptional numerical seeds
   are separate features; the vertex change must not require them all at once.

Use small cases in existing relevant test files. Inspect type stability and
measure allocation/runtime changes against the current tableau traversal.
Preserve opportunities for independent per-weight numerical work without
rewriting the existing traversal solely to parallelize it. Apply the
[shared implementation policy](crystal_basis_integration_plan.md#shared-implementation-policy).

Completion means the existing construction can use tensor-word vertices in place
of tableaux, fundamental embeddings compose correctly, and consumers use the
same weight/raising/lowering operations. Part 2 consumes these vertex operations
directly; it does not require owned whole-crystal graphs or precomputed edge tables.

## References and provenance


Primary reference: Daniel Bump and Anne Schilling, *Crystal Bases:
Representations and Combinatorics* (2017), the supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

- §2.2, Definition 2.13 and equation (2.6): crystal axioms and string lengths.
- §§2.3–2.4, printed pp. 18–23: tensor convention, associativity, and signatures.
- §5.7, Theorem 5.20 and its proof, printed pp. 72–73: normal crystals and the
  highest component in tensor products of fundamental crystals.
- Chapters 3 and 6: tableau realizations retained as seed builders/references.

These passages were consulted in the local PDF's extracted text for this rewrite.
The vertex type and storage choices are engineering proposals, not prescriptions from
the book. This documentation-only change ran no numerical tests or doc builds.
