# Integer q-labels and weight coordinates for E6, E7, E8, and F4

Status: coordinate specification. E6/E7 conversions are implemented in
[E6.jl](../src/Symmetry/E6.jl) and [E7.jl](../src/Symmetry/E7.jl);
[F4.jl](../src/Symmetry/F4.jl) also implements the F4 conversion below.
E8 remains a coordinate proposal.
Existing A–D/G2 weight-coordinate conventions remain unchanged.

## Storage and ordering

Keep q-labels in Dynkin-node order, but store **reversed scaled simple-root
coordinates**. For a weight ``\mu=\sum_i c_i\alpha_i``, define

```math
z=dPc,\qquad (Pc)_j=c_{r+1-j},\qquad P^2=I.
```

The reversal applies only to stored weight coordinates, not to q-labels,
Cartan nodes, or lowering-operator indices. All vectors below are columns;
Julia storage remains `NTuple{r,Int}`.

| Family | Rank r | Scale d | Stored coordinate |
| --- | ---: | ---: | --- |
| E6 | 6 | 3 | ``z_j=3c_{7-j}`` |
| E7 | 7 | 2 | ``z_j=2c_{8-j}`` |
| E8 | 8 | 1 | ``z_j=c_{9-j}`` |
| F4 | 4 | 1 | ``z_j=c_{5-j}`` |

The previous unreversed proposal made every lowering decrease the reversed
lexicographic order, but did not match the required operator-to-slot order.
With this convention, operator 1 reduces the last coordinate; operator 2
leaves the last fixed and reduces the second-last; operator i reduces slot
r+1−i and leaves every other coordinate fixed. Thus every nonzero lowering
strictly decreases `reverse(z)` lexicographically, as required by `rev_less`
and `sorted_zvals` in [sort.jl](../src/sort.jl).

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

## Conversion tables

With the Cartan convention above, Dynkin coordinates satisfy ``a=Ac``.
Define ``M=dPA^{-1}``. Then

```math
\operatorname{qlab2mwz}(q)=Mq,\qquad
\operatorname{getqlabel}(z)=APz/d,\qquad
\operatorname{getdz}(i)=d\mathbf e_{r+1-i}.
```

The factors d are the least denominators clearing the inverse Cartan matrices.
Column k of each M below is the maximal stored weight of fundamental k.
These are row reversals of the previous tables; do not reverse their columns.

### E6

```math
M_{E6}=\begin{pmatrix}
2 & 3 & 4 & 6 & 5 & 4 \\
4 & 6 & 8 & 12 & 10 & 5 \\
6 & 9 & 12 & 18 & 12 & 6 \\
5 & 6 & 10 & 12 & 8 & 4 \\
3 & 6 & 6 & 9 & 6 & 3 \\
4 & 3 & 5 & 6 & 4 & 2
\end{pmatrix}.
```

### E7

```math
M_{E7}=\begin{pmatrix}
2 & 3 & 4 & 6 & 5 & 4 & 3 \\
4 & 6 & 8 & 12 & 10 & 8 & 4 \\
6 & 9 & 12 & 18 & 15 & 10 & 5 \\
8 & 12 & 16 & 24 & 18 & 12 & 6 \\
6 & 8 & 12 & 16 & 12 & 8 & 4 \\
4 & 7 & 8 & 12 & 9 & 6 & 3 \\
4 & 4 & 6 & 8 & 6 & 4 & 2
\end{pmatrix}.
```

### E8

```math
M_{E8}=\begin{pmatrix}
2 & 3 & 4 & 6 & 5 & 4 & 3 & 2 \\
4 & 6 & 8 & 12 & 10 & 8 & 6 & 3 \\
6 & 9 & 12 & 18 & 15 & 12 & 8 & 4 \\
8 & 12 & 16 & 24 & 20 & 15 & 10 & 5 \\
10 & 15 & 20 & 30 & 24 & 18 & 12 & 6 \\
7 & 10 & 14 & 20 & 16 & 12 & 8 & 4 \\
5 & 8 & 10 & 15 & 12 & 9 & 6 & 3 \\
4 & 5 & 7 & 10 & 8 & 6 & 4 & 2
\end{pmatrix}.
```

### F4

```math
M_{F4}=\begin{pmatrix}
2 & 4 & 3 & 2 \\
4 & 8 & 6 & 3 \\
3 & 6 & 4 & 2 \\
2 & 3 & 2 & 1
\end{pmatrix}.
```

## Inverse, lattice, and lowering checks

Let ``u=Pz`` (equivalently, `u = reverse(z)`). For E6–E8,

```math
a_j=(2u_j-\sum_{k\sim j}u_k)/d.
```

For F4,

```math
\begin{aligned}
a_1&=2z_4-z_3,\\
a_2&=-z_4+2z_3-z_2,\\
a_3&=-2z_3+2z_2-z_1,\\
a_4&=-z_2+2z_1.
\end{aligned}
```

For E6/E7, require exact divisibility of every entry of `APz` by d;
an arbitrary integer tuple need not represent an integral weight. E8/F4
accept every integer tuple as a lattice weight. Dominance requires
``APz/d\geq0`` componentwise, not merely nonnegative z. The inverse gives
an irrep q-label only for a highest weight; for other weights it gives signed
Dynkin coordinates.

`getdz` is the quantity subtracted, not the signed increment. Lowering and
raising at node i subtract and add ``d\mathbf e_{r+1-i}``, respectively.
In Dynkin coordinates, lowering subtracts column i of A. The identities

```math
APM=dI,\qquad MA=dP,\qquad
M(a-A_{:,i})=Ma-d\mathbf e_{r+1-i}
```

make these descriptions consistent. The crystal pairing is
``\phi_i-\epsilon_i=(APz)_i/d``. These formulas specify weight changes;
the crystal or numerical representation determines whether an action exists.

| Family and highest q-label | Maximal z | First lowering | Resulting z |
| --- | --- | --- | --- |
| E6, `(1,0,0,0,0,0)` | `(2,4,6,5,3,4)` | color 1 | `(2,4,6,5,3,1)` |
| E7, `(0,0,0,0,0,0,1)` | `(3,4,5,6,4,3,2)` | color 7 | `(1,4,5,6,4,3,2)` |
| E8, `(0,0,0,0,0,0,0,1)` | `(2,3,4,5,6,4,3,2)` | color 8 | `(1,3,4,5,6,4,3,2)` |
| F4, `(0,0,0,1)` | `(2,3,2,1)` | color 4 | `(1,3,2,1)` |

The final F4 tuple has Dynkin coordinates `(0,0,1,-1)`, as expected for a
non-highest weight. An example's first available lowering need not be color 1.

## Cartan operators and root-space metric

The stored coordinates are not orthogonal Cartesian coordinates. They can be
interpreted as eigenvalues of

```math
Z_j=d\sum_k(PA^{-1})_{jk}H_k,\qquad
[Z_j,F_i]=-d\delta_{j,r+1-i}F_i,
```

where ``H_k=\alpha_k^\vee``. Recover coroot eigenvalues as ``APz/d``.
For scalar products, let ``D_{ii}=(\alpha_i,\alpha_i)/2`` and ``G=DA``.
Then

```math
(\mu,\nu)=z_\mu^{\mathsf T}PGPz_\nu/d^2.
```

Use D=I for E6–E8 and ``D=\operatorname{diag}(1,1,1/2,1/2)`` for F4.
This coordinate metric is distinct from the representation-basis Gram matrices
in `Irep.innerprod`. Weight addition and negation remain componentwise.

## Implementation guidance and validation

Precompute M and the root-step tuples. Use short integer sums and exact division
for the inverse; no matrix inversion or dense intermediates belong in per-weight
loops. Keep tuple returns concrete and check arithmetic overflow. No new type
is needed. Do not store Dynkin coordinates directly as z: lowering can increase
a later coordinate and violate the existing sorting order.

Validate with exact arithmetic:

- `APM = dI`, `MA = dP`, and integral M with the stated minimal scale.
- Fundamental and mixed-label conversion round trips.
- Every color's lowering identity, unchanged slots, and positive root-step
  entry at precisely slot r+1−i.
- Strict decrease under reversed-coordinate lexicographic comparison.
- Symmetry of `PGP` and agreement with the unreversed root-space metric.

These are coordinate checks, not exceptional-irrep or CG tests. The linked
references specify roots and node numbering; the storage reversal is a LurCGT
convention. No Julia implementation or documentation build is required for this
proposal change.
