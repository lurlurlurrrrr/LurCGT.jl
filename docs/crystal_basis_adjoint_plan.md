# Part 3c: E8 adjoint crystal and numerical seed

Status: design proposal. Uses the recipe/extraction machinery of
[part 3a](crystal_basis_minuscule_plan.md) and the
[shared requirements](crystal_basis_integration_plan.md).
This work can proceed independently of part 3b once that machinery is ready.

## Scope and sequence

1. Enumerate E8 roots with the shared Cartan/node convention.
2. Construct the adjoint crystal of §5.5: one vertex per root and eight distinct
   zero-weight vertices. Define reversible crystal operators and string lengths.
3. Construct the numerical adjoint action in a Chevalley basis, including
   structure constants and the zero-weight metric. Crystal adjacency is not a
   substitute for these matrices.
4. Connect the 248-dimensional seed to the incremental fundamental recipes below
   and part 3a's general-highest-weight construction.

## Fundamental dependencies

E8 has no minuscule seed. Set A=``B(\omega_8)`` (the 248-dimensional adjoint).
The book gives explicit highest tuples, equivalently the following incremental
recipes. Here ``v_\beta`` is the adjoint vertex of root weight β, and ``u_j``
is the selected highest vertex of ``B(\omega_j)``. Use the book's tensor order.

| Target node | Highest tuple generating its component | Dimension |
| --- | --- | ---: |
| 8 | Adjoint highest root vertex ``v_{\omega_8}`` | 248 |
| 1 | ``v_{\omega_1-\omega_8}\otimes u_8`` | 3,875 |
| 7 | ``v_{\omega_7-\omega_8}\otimes u_8`` | 30,380 |
| 2 | ``v_{\omega_2-\omega_1}\otimes u_1`` | 147,250 |
| 6 | ``v_{\omega_6-\omega_7}\otimes u_7`` | 2,450,240 |
| 3 | ``v_{\omega_3-\omega_2}\otimes u_2`` | 6,696,000 |
| 5 | ``v_{\omega_5-\omega_6}\otimes u_6`` | 146,325,270 |
| 4 | ``v_{\omega_4-\omega_5}\otimes u_5`` | 6,899,079,264 |

Each displayed difference is a root, with the string lengths listed in Table
5.4. Verify that all raising operators kill each tuple before enumeration.
This records a construction for every node without claiming eager storage of
all E8 fundamentals is practical. Node 4 alone has billions of vertices.

## Completion

Validate 240 root vertices and eight zero-weight vertices, the highest-root
label, crystal inverses, numerical Chevalley/Serre relations, and the metric.
Check each recipe's highest tuple before enumeration. Exercise manageable
components and keep all recipes queryable without allocating their full spaces.
Complete recipe coverage does not imply practical eager construction of the
largest E8 representations.

## Reference

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
user-supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

Relevant passages: §5.5 and §5.6, Table 5.4.
