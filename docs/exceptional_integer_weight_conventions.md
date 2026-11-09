# Integer q-labels and weight coordinates for E6, E7, E8, and F4

Status: standalone coordinate proposal. This adds no symmetry implementation,
representation backend, or spin support, and does not change existing A–D/G2
conventions or the scope of the crystal-vertex replacement.

## Recommendation

Yes: q-labels and every weight can be represented by rank-length integer tuples
for all four families. Keep q-labels as fundamental-weight coefficients (Dynkin
labels), and store weights in **scaled simple-root coordinates**. This gives
explicit integral conversions and particularly simple lowering steps while
preserving LurCGT's current highest-first sorting convention.

Let q be the q-label of the highest weight λ, A the Cartan matrix defined below,
and d the family-dependent denominator-clearing factor. Define ``M=dA^{-1}``.
Then the proposed API meanings are

```math
\boxed{\operatorname{qlab2mwz}(q)=Mq,\qquad
\operatorname{getqlabel}(z)=Az/d,\qquad
\operatorname{getdz}(i)=d\,\mathbf e_i.}
```

All vectors in the formulas are columns; Julia can continue storing weights and
q-labels as `NTuple{r,Int}`. Here ``\mathbf e_i`` is the coordinate unit vector,
not a crystal raising operator.

| Family | Rank r | d | Stored weight z for ``\mu=\sum_j c_j\alpha_j`` |
| --- | ---: | ---: | --- |
| E6 | 6 | 3 | ``z_j=3c_j`` |
| E7 | 7 | 2 | ``z_j=2c_j`` |
| E8 | 8 | 1 | ``z_j=c_j`` |
| F4 | 4 | 1 | ``z_j=c_j`` |

These factors are the least common denominators of the corresponding inverse
Cartan matrices, computed exactly below. They apply to **every integral weight**,
not only dominant highest weights. The stored coordinates are not orthogonal
Cartesian coordinates; their metric is given later in this note.

## Node numbering and Cartan convention

Define

```math
A_{ij}=\langle\alpha_i^\vee,\alpha_j\rangle,
\qquad
\langle\alpha_i^\vee,\omega_j\rangle=\delta_{ij}.
```

Consequently the Dynkin-label column of ``\alpha_i`` is column i of A.
This fixes the transpose convention, which matters for F4. This convention is
also stated explicitly in [Sage's Lie group basics](https://doc.sagemath.org/html/en/thematic_tutorials/lie/lie_basics.html#the-cartan-matrix).

For E6–E8, use the chain `1–3–4–5–6[–7[–8]]` with node `2` attached to `4`.
Thus A has diagonal 2 and entries −1 on each edge, zero elsewhere. These are
[Sage's documented E-family node labels](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_E.html#sage.combinat.root_system.type_E.CartanType.dynkin_diagram).
The common edge set is `(1,3), (2,4), (3,4), (4,5), (5,6)`; append `(6,7)`
for E7 and also `(7,8)` for E8.

For F4, nodes 1,2 are long and nodes 3,4 short, with chain `1–2⇒3–4`. The
simple roots in orthonormal coordinates can be chosen as

```math
\alpha_1=e_2-e_3,\quad \alpha_2=e_3-e_4,\quad
\alpha_3=e_4,\quad \alpha_4=(e_1-e_2-e_3-e_4)/2.
```

These are the [documented F4 simple roots](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_F.html#sage.combinat.root_system.type_F.AmbientSpace.simple_root).
With the pairing convention above they give

```math
A_{F4}=\begin{pmatrix}
2&-1&0&0\\
-1&2&-1&0\\
0&-2&2&-1\\
0&0&-1&2
\end{pmatrix}.
```

In particular ``A_{23}=-1`` and ``A_{32}=-2``. Do not copy a transposed
Cartan matrix while retaining the same conversion formulas.

## Q-label to maximal weight, and the inverse

For the irreducible highest weight, take

```math
q=(q_1,\ldots,q_r),\quad q_i\in\mathbb Z_{\geq0},\qquad
\lambda=\sum_i q_i\omega_i.
```

This labels Lie-algebra irreps, equivalently irreps of the simply connected
compact group; a specified quotient group can restrict which labels are allowed.
The conventions below do not impose such a quotient restriction.

Write any weight as ``\mu=\sum_i c_i\alpha_i`` and let
``a_i=\langle\alpha_i^\vee,\mu\rangle``. Then ``a=Ac``. Set ``z=dc``,
so ``z=Ma`` and ``a=Az/d``. For the highest weight, ``a=q``. This proves
both conversion directions and explains why the same conversion works on
all weights, even though their Dynkin coordinates need not be nonnegative.

The following matrices are derived by exact inversion of the Cartan matrices
specified above. Column k of M is the stored maximal weight for q-label
``q=\mathbf e_k``. They are proposed LurCGT coordinate tables, not copied
orthogonal-coordinate weight tables.

### E6

```math
M_{E6}=\begin{pmatrix}
4 & 3 & 5 & 6 & 4 & 2 \\
3 & 6 & 6 & 9 & 6 & 3 \\
5 & 6 & 10 & 12 & 8 & 4 \\
6 & 9 & 12 & 18 & 12 & 6 \\
4 & 6 & 8 & 12 & 10 & 5 \\
2 & 3 & 4 & 6 & 5 & 4
\end{pmatrix}.
```

### E7

```math
M_{E7}=\begin{pmatrix}
4 & 4 & 6 & 8 & 6 & 4 & 2 \\
4 & 7 & 8 & 12 & 9 & 6 & 3 \\
6 & 8 & 12 & 16 & 12 & 8 & 4 \\
8 & 12 & 16 & 24 & 18 & 12 & 6 \\
6 & 9 & 12 & 18 & 15 & 10 & 5 \\
4 & 6 & 8 & 12 & 10 & 8 & 4 \\
2 & 3 & 4 & 6 & 5 & 4 & 3
\end{pmatrix}.
```

### E8

```math
M_{E8}=\begin{pmatrix}
4 & 5 & 7 & 10 & 8 & 6 & 4 & 2 \\
5 & 8 & 10 & 15 & 12 & 9 & 6 & 3 \\
7 & 10 & 14 & 20 & 16 & 12 & 8 & 4 \\
10 & 15 & 20 & 30 & 24 & 18 & 12 & 6 \\
8 & 12 & 16 & 24 & 20 & 15 & 10 & 5 \\
6 & 9 & 12 & 18 & 15 & 12 & 8 & 4 \\
4 & 6 & 8 & 12 & 10 & 8 & 6 & 3 \\
2 & 3 & 4 & 6 & 5 & 4 & 3 & 2
\end{pmatrix}.
```

### F4

```math
M_{F4}=\begin{pmatrix}
2 & 3 & 2 & 1 \\
3 & 6 & 4 & 2 \\
4 & 8 & 6 & 3 \\
2 & 4 & 3 & 2
\end{pmatrix}.
```

The E-family inverse can be written without a dense matrix multiply:

```math
q_j=\frac{2z_j-\sum_{k\sim j}z_k}{d},
```

where the sum runs over neighbors in the specified Dynkin diagram. For F4,
where d=1, the inverse is explicitly

```math
\begin{aligned}
q_1&=2z_1-z_2,\\
q_2&=-z_1+2z_2-z_3,\\
q_3&=-2z_2+2z_3-z_4,\\
q_4&=-z_3+2z_4.
\end{aligned}
```

For E6/E7, an arbitrary integer tuple z is a valid integral weight **only if**
every entry of Az is divisible by d. Equivalently, the stored weight lattice
is ``M\mathbb Z^r``. Use exact divisibility checks, not truncating division.
For E8/F4, every integer z is a lattice weight. In all cases, dominance requires
``Az/d\geq0`` componentwise; nonnegative z alone does not imply dominance.

`getqlabel` on a known maximal weight returns its irrep q-label. Applying the
same algebraic inverse to a non-highest weight gives signed Dynkin coordinates,
not the q-label of a different irrep. Keep that distinction in callers.

## Decrease under each lowering operator

Whenever a lowering operator acts nontrivially, its weight changes by
``\mu\mapsto\mu-\alpha_i``. Therefore

```math
z\mapsto z-d\mathbf e_i.
```

| Family | `getdz(i)` | Actual change to stored weight |
| --- | --- | --- |
| E6 | 3 in slot i, zero elsewhere | ``z_i\gets z_i-3`` |
| E7 | 2 in slot i, zero elsewhere | ``z_i\gets z_i-2`` |
| E8 | 1 in slot i, zero elsewhere | ``z_i\gets z_i-1`` |
| F4 | 1 in slot i, zero elsewhere | ``z_i\gets z_i-1`` |

This uses the same sign convention as the current code: `getdz` returns the
quantity **subtracted**, as in `nw = ow .- getdz(S, i)` in
[`src/irep.jl`](../src/irep.jl). It is not the signed increment `nw - ow`.
Raising adds the same tuple. After k steps of color i, subtract ``kd\mathbf e_i``.

In Dynkin coordinates the same operation is ``a\mapsto a-A_{:,i}``;
``M A=dI`` makes these descriptions identical. These formulas determine the
weight change, not whether an operator is nonzero on a particular vertex or
vector. Existence of the action still comes from the crystal or numerical module.
The crystal identity is ``\phi_i-\epsilon_i=(Az)_i/d``.

### Examples

| Family and highest q-label | Maximal z | First lowering | Resulting z |
| --- | --- | --- | --- |
| E6, `(1,0,0,0,0,0)` | `(4,3,5,6,4,2)` | color 1 | `(1,3,5,6,4,2)` |
| E7, `(0,0,0,0,0,0,1)` | `(2,3,4,6,5,4,3)` | color 7 | `(2,3,4,6,5,4,1)` |
| E8, `(0,0,0,0,0,0,0,1)` | `(2,3,4,6,5,4,3,2)` | color 8 | `(2,3,4,6,5,4,3,1)` |
| F4, `(0,0,0,1)` | `(1,2,3,2)` | color 4 | `(1,2,3,1)` |

For the F4 example, the last tuple converts to Dynkin coordinates `(0,0,1,-1)`.
The negative entry is expected for a non-highest weight.

## Fit with LurCGT's existing conventions

The current library already distinguishes q-labels from stored z-weights:
[`SU.jl`](../src/Symmetry/SU.jl), [`Sp.jl`](../src/Symmetry/Sp.jl),
[`SO.jl`](../src/Symmetry/SO.jl), and [`G2.jl`](../src/Symmetry/G2.jl) supply
family-specific `qlab2mwz`, `getqlabel`, and `getdz`. For example, G2 maps
`(a,b)` to `(a,a+2b)`; SO uses doubled orthogonal coordinates. Thus exceptional
families can choose a different integral coordinate basis without changing the
meaning of the existing families' tuples.

This proposal keeps:

- `nlops = nzops = rank`, with rank-length integer tuple sector keys.
- Addition of weights in tensor products and negation under dualization.
- A constant `getdz(i)` for each simple lowering operator.
- The highest-first reverse-lexicographic sort in
  [`src/sort.jl`](../src/sort.jl), also used by
  [`src/clebsch.jl`](../src/clebsch.jl). Every lowering strictly decreases one
  coordinate and leaves the others fixed, so it strictly decreases that order.

One could instead store Dynkin coordinates directly, making `qlab2mwz(q)=q` and
`getdz(i)=A[:,i]`. However, off-diagonal negative entries cause some coordinates
to increase under lowering. For example, in F4 subtracting column 1 increases
the second coordinate. The existing reverse-lexicographic sort would then place
that lowered weight above its source. Scaled root coordinates avoid this ordering
change without introducing an extra grade or extra stored coordinate.

The z values can also be interpreted as eigenvalues of commuting Cartan operators
``Z_j=d\sum_k(A^{-1})_{jk}H_k``, where ``H_k=\alpha_k^\vee``. Then
``[Z_j,F_i]=-d\delta_{ji}F_i``. For Chevalley relations, recover the usual coroot
eigenvalues as ``H_i:(Az)_i/d``; z itself is not the list of coroot eigenvalues.

If a root-space scalar product is needed, use the coordinate metric rather than
an ordinary dot product. Let ``D_{ii}=(\alpha_i,\alpha_i)/2``, so ``G=DA`` is
the symmetric simple-root Gram matrix. Then

```math
(\mu,\nu)=\frac{1}{d^2}z_\mu^{\mathsf T}Gz_\nu.
```

For E6–E8, take D=I; for F4's stated normalization, take
``D=\operatorname{diag}(1,1,1/2,1/2)``. This root-space metric is distinct from
the existing Gram matrices of representation basis vectors in `Irep.innerprod`.

## Implementation guidance and checks

Precompute the small integer M matrices and fixed root-step tuples. Do not invert
A or allocate dense matrix intermediates in per-weight loops. The inverse
conversion uses short integer sums and exact division. Keep tuple return types
concrete; use checked arithmetic or exact intermediates with checked conversion
when labels can overflow `Int`. No new data type is needed for these conversions.
Defining new exceptional symmetry types or numerical backends remains separate work.

The following checks were run for this note with exact rational/integer arithmetic:

- ``AM=MA=dI`` and the denominator-clearing factors for all four families.
- Conversion round trips for every fundamental q-label and one mixed label per
  family; the lowering identity ``M(a-A_{:,i})=Ma-d\mathbf e_i`` for every color.
- Strict decrease under the existing reversed-coordinate lexicographic order.
- F4's Cartan matrix independently reconstructed from the stated Euclidean roots.
- Weyl-reflection closure of the simple roots: 72, 126, 240, and 48 roots for
  E6, E7, E8, and F4 respectively, with each root having uniformly nonnegative
  or nonpositive simple-root coordinates.

These are coordinate/algebra checks, not tests of exceptional representations or
CG coefficients. No Julia implementation, DMRG test, or documentation build was run.
The coordinate choice and matrices here are derived recommendations; the linked
Sage documentation supplies the root conventions and node numbering.
