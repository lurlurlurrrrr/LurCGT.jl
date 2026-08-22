# CGT and related symbols

A Clebsch-Gordan tensor (CGT) is an intertwiner between products of irreducible
representations (irreps). LurCGT uses CGTs to separate symmetry-determined
coefficient data from the reduced tensors that carry a tensor-network
calculation's variational data.

## Arbitrary-rank CGTs

A CGT may have any number of input and output legs. Its physical axes carry
irrep basis states, and each leg has a q-label and direction. The familiar
three-leg coupling is the special case
``q_1 \otimes q_2 \to q_3``. More generally, a CGT represents an invariant map
between two products of irreps, such as
``q_1 \otimes \cdots \otimes q_m \to r_1 \otimes \cdots \otimes r_n``.

LurCGT stores a CGT block-sparsely by weight sector. The final array axis is an
outer-multiplicity axis; it is separate from the physical irrep-basis axes.
For an ordinary three-leg coupling, use `getNsave_cg3`. Internally, higher-rank
CGTs are assembled and transformed through fusion trees and the cached
recoupling data described below.

## Outer multiplicity and the canonical basis

The same product of irreps can contain an output irrep more than once. The
number of independent intertwiners is its outer multiplicity. A CGT therefore
requires an additional label to distinguish these independent coupling
channels, even after all external q-labels have been fixed.

LurCGT fixes this freedom by using a canonical fusion-tree basis. The basis is
determined by the ordered external q-labels, intermediate irreps, and
outer-multiplicity labels of the constituent three-leg couplings. The
bookkeeping object returned by `get_CGTom` describes the flattened canonical
outer-multiplicity space. This fixed convention makes cached CGTs and basis
transforms reproducible.

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
  are contracted. `getNsave_Xsymbol` returns the dense coefficient array for
  the requested contracted legs.

These symbols are used internally to keep contractions, permutations, and
tensor decompositions in the canonical CGT basis. They are cached because the
underlying coefficients are exact but can be expensive to construct.

## Related CGT transformations

`getNsave_CGTSVD` and `getNsave_CGTQR` construct the CGT-side basis changes
needed by symmetry-respecting SVD and QR decompositions. They leave the
reduced-tensor calculation to the caller while accounting for the symmetry and
outer-multiplicity structure of the split. `to_float` converts an exact CGT to
a numerical sparse array when an explicit coefficient representation is needed.
