# LurCGT.jl

LurCGT provides exact Clebsch-Gordan and representation-theory data for
non-Abelian symmetries. It is the symmetry backend used by
[Telum.jl](https://github.com/ssblee/Telum.jl). It may later be used by other
tensor libraries.

The package computes representation data deterministically with exact arithmetic
and stores generated results in an SQLite cache. Stored data can be reused
arbitrarily.

```@contents
Pages = ["getting-started.md", "cgt-and-related-symbols.md", "api.md"]
Depth = 2
```
