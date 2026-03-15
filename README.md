# LurCGT.jl

`LurCGT.jl` is the Clebsch-Gordan and symmetry core extracted from the original
`CGTfromInts` repository. It provides symmetry types, irreps, CGT/F/R/X-symbol
machinery, decomposition helpers, and the SQLite-backed storage layer.

## Installation

Once this repository is pushed to its own remote, install it with:

```julia
using Pkg
Pkg.add(url="<LurCGT-repo-url>")
```

For local development:

```julia
using Pkg
Pkg.develop(path=".")
```

## Split Notes

`QSpace` functionality now lives in `QSpaces.jl`. Typical post-split usage is:

```julia
using LurCGT
using QSpaces
```

## Testing

```julia
using Pkg
Pkg.test()
```
