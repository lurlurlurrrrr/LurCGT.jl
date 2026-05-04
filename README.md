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

## SQLite Database Directories

By default, `LurCGT.jl` stores SQLite data under `~/.LurCGT_sqlite`.

- Global databases default to `~/.LurCGT_sqlite/global`
- Local per-process databases default to `~/.LurCGT_sqlite/local`
- Merge lock files default to `~/.LurCGT_sqlite/global/locks`

You can override these locations with environment variables before using the package:

- `LURCGT_RUN_MODE`: selects which directory set to use; supported values are `local` and `server`
- `LURCGT_LOCALDB_DIR`: local-mode directory for per-process local databases
- `LURCGT_GLOBALDB_DIR`: local-mode directory for global `*.db` files
- `LURCGT_LOCALDB_DIR_NODE`: server-mode directory for per-process local databases
- `LURCGT_GLOBALDB_DIR_NODE`: server-mode directory for global `*.db` files

`LURCGT_RUN_MODE` defaults to `local`.

In `server` mode, `LurCGT.jl` uses the `_NODE` directories. If a node DB does not
exist yet, it first looks for the matching DB in the non-`_NODE` directory and copies
it into the `_NODE` directory. If no non-`_NODE` DB exists, it creates a fresh DB in
the `_NODE` directory.

## Testing

```julia
using Pkg
Pkg.test()
```

