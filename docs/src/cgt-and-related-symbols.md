# CGT and related symbols

A Clebsch-Gordan tensor (CGT) generalizes ordinary three-dimensional
Clebsch-Gordan coefficients. LurCGT stores and uses a CGT as a tensor of
symmetry-determined coefficients: its ordinary axes are basis-state axes for the
participating irreducible representations (irreps), and its last axis
enumerates outer-multiplicity channels. This coefficient tensor separates
symmetry data from the reduced tensors that carry a tensor-network
calculation's variational data.

Mathematically, the same data can also be described as an intertwiner between
products of irreps. This page emphasizes the tensor viewpoint because it is the
form used by LurCGT's storage, contraction, and basis-change routines.

## Arbitrary-rank CGTs

A CGT may have any number of input and output legs. Its physical axes carry
irrep basis states, and each leg has a q-label and direction. The familiar
three-leg coupling is the special case
``q_1 \otimes q_2 \to q_3``. More generally, a CGT has one tensor leg for each
input or output irrep, such as
``q_1, \ldots, q_m, r_1, \ldots, r_n`` for the coupling conventionally written
as ``q_1 \otimes \cdots \otimes q_m \to r_1 \otimes \cdots \otimes r_n``.
LurCGT's stored CGT convention orders input legs first and output legs second;
within each group, legs are sorted by ascending q-label.

LurCGT stores the basic three-leg Clebsch-Gordan coefficient tensors and uses
them as building blocks for higher-rank CGTs. Internally, arbitrary-rank CGTs are
assembled and transformed through fusion trees from the stored three-leg
coefficients and the cached F- and R-symbol data described below; they are not
stored as independent primitive CGT objects.

## Outer multiplicity and the canonical basis

The same product of irreps can contain an output irrep more than once. The
number of independent intertwiners is its outer multiplicity. A CGT therefore
requires an additional label to distinguish these independent coupling
channels, even after all external q-labels have been fixed.

LurCGT fixes this freedom by using a canonical fusion-tree basis. The basis is
determined by the ordered external q-labels, intermediate irreps, and
outer-multiplicity labels of the constituent three-leg couplings. This fixed
convention makes cached CGTs and basis transforms reproducible.

Changing the fusion tree or permuting legs can mix outer-multiplicity channels.
Use `getNsave_CGTperm` for the corresponding canonical-basis transformation;
`get_conj_perm` gives the related transformation for conjugation.

## F-, R-, and X-symbols

The recoupling symbols express elementary changes of CGT basis.

- **F-symbols** implement associativity: they change between the two fusion
  trees for ``(q_1 \otimes q_2) \otimes q_3`` and
  ``q_1 \otimes (q_2 \otimes q_3)``. Obtain the cached exact matrix with
  `getNsave_Fsymbol`.
- **R-symbols** implement exchange for the supported equal-input fusion
  ``q \otimes q \to r``. Obtain the cached exact matrix with
  `getNsave_Rsymbol`.
- **X-symbols** recouple the outer-multiplicity bases that result when two CGTs
  are contracted. They are computed from sequences of F- and R-symbol moves.
  `getNsave_Xsymbol` returns the dense coefficient array for the requested
  contracted legs.

These symbols are used internally to keep contractions, permutations, and
tensor decompositions in the canonical CGT basis. The primitive stored data are
the basic three-leg Clebsch-Gordan coefficients together with F- and R-symbols.

## Related CGT transformations

`getNsave_CGTSVD` and `getNsave_CGTQR` construct the CGT-side basis changes
needed by symmetry-respecting SVD and QR decompositions. They leave the
reduced-tensor calculation to the caller while accounting for the symmetry and
outer-multiplicity structure of the split. `to_float` converts an exact CGT to
a numerical sparse array when an explicit coefficient representation is needed.

## Further reading

This page gives only the conventions needed to use LurCGT. For more background,
see:

- Clebsch-Gordan tensor:
  [SciPost Physics Codebases 40](https://scipost.org/SciPostPhysCodeb.40).
- X-symbols for non-Abelian tensor contractions:
  [Phys. Rev. Research 2, 023385](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023385).
- F- and R-symbols in tensor-category notation:
  [TensorKit.jl category appendix](https://quantumkithub.github.io/TensorKit.jl/stable/appendix/categories/).
