# Part 3b: Irreps for new symmetries using virtual crystals

Status: design proposal. Depends on [part 3a](crystal_basis_minuscule_plan.md)
and follows the [shared requirements](crystal_basis_integration_plan.md).

## Scope and sequence

1. Use the E6 minuscule seeds from part 3a to construct the 26-vertex F4 seed
   through the book's virtual operators.
2. Supply its exact numerical `Irep` realization; the virtual crystal alone
   does not define classical lowering matrices.
3. Integrate the seed with part 3a's recipes for all F4 fundamentals and general
   highest weights, including weight adapters, duality, and cache storage.
4. Use the virtual D4 realization of G2 as a comparison with the existing G2
   backend. Preserve that numerical seed. D4 spin factors require part 3a.

## Virtual operators

For virtualization record ambient type, node orbits ``O_i``, scale factors
``\gamma_i``, and

```math
\Psi(\omega_i)=\gamma_i\sum_{j\in O_i}\widehat\omega_j,
\qquad
\widetilde f_i^{\rm virt}=\prod_{j\in O_i}
\widehat{\widetilde f}_j^{\gamma_i}.
```

Use the analogous raising operators. Orbit nodes commute. At every retained
vertex check equality of ambient string lengths within each orbit, divisibility
by ``\gamma_i``, closure under both virtual operators, and the weight map.
These are the alignment and normality requirements in §§5.1–5.3, not optional
checks on a visually folded graph.

The book's explicit exceptional constructions (§5.6, printed p. 70) are:

| Target seed | Ambient crystal and highest vertex | Virtual lowering operators in target node order |
| --- | --- | --- |
| F4, ``B(\omega_4)`` (26) | E6 ``B(\widehat\omega_1)\otimes B(\widehat\omega_6)``, highest⊗highest | ``\widehat f_2^2,\ \widehat f_4^2,\ \widehat f_3\widehat f_5,\ \widehat f_1\widehat f_6`` |
| G2, ``B(\omega_1)`` (7) | D4 ``B(\widehat\omega_1)\otimes B(\widehat\omega_3)\otimes B(\widehat\omega_4)``, all highest | ``\widehat f_1\widehat f_3\widehat f_4,\ \widehat f_2^3`` |

Hats in this table mean **crystal** operators. In F4 the factors are
``(2,2,1,1)``; in G2 they are ``(1,3)``. These are the book's similarity
embeddings, not naive restriction/branching of the 27 or 8 dimensional module.
In particular, a product of classical ambient lowering matrices is generally
not a target Chevalley generator: its weight change and commutators differ.

The separate F4 numerical realization and its acceptance gate are specified in
[numerical realization](#numerical-realization).
Keep the existing validated G2 numerical seed.

## Numerical realization

For F4 numerical matrices, a separate explicit realization is therefore required.
A general exact highest-weight-module construction is a suitable fallback:
use Chevalley–Serre relations, a PBW ordering and the highest-weight relations
``E_i v_\lambda=0``, ``H_i v_\lambda=\lambda_i v_\lambda``,
``F_i^{\lambda_i+1}v_\lambda=0``; compute the contravariant form and quotient
its radical weight by weight. This is substantial work and must be validated
before claiming virtual crystals alone supply a new `Irep` backend. A specialized
26-dimensional F4 realization can replace this fallback if its normalization and
relations are proved. Keep the existing validated G2 numerical seed.

## Completion

Check virtual-operator inverses, alignment, weight conversion, and the expected
F4/G2 seed dimensions. Validate the F4 numerical generators and metric under
the target Cartan convention, then construct the other F4 fundamentals through
the shared recipes. Run existing CG/recoupling checks for enabled numerical
representations. Do not treat ambient classical operator products as the
numerical implementation of virtual crystal operators.

## Reference

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
user-supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

Relevant passages: §§5.1–5.3 and §5.6, printed p. 70.
