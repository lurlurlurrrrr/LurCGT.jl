# LurCGT.jl

LurCGT provides exact Clebsch-Gordan and representation-theory data for
non-Abelian symmetries. It is the symmetry backend used by
[Telum.jl](https://github.com/ssblee/Telum.jl), and can also be used directly
to construct symmetry-adapted tensors.

The package computes representation data with exact arithmetic and stores
generated results in an SQLite cache. This makes coefficient generation
reproducible and avoids accumulating floating-point error while constructing
CGT, F-, R-, and X-symbol data.

Most tensor-network applications should use LurCGT through Telum. The public
API documented here is stable for direct integrations.

```@contents
Pages = ["getting-started.md", "api.md"]
Depth = 2
```
