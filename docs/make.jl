using Documenter
using LurCGT

makedocs(
    modules = [LurCGT],
    sitename = "LurCGT.jl",
    authors = "Kiyeon Kim",
    repo = "github.com/lurlurlurrrrr/LurCGT.jl/blob/{commit}{path}#{line}",
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "CGT and related symbols" => "cgt-and-related-symbols.md",
        "Public API" => [
            "Overview" => "api.md",
            "Symmetries" => "api/symmetries.md",
            "Fusion and CGT transforms" => "api/cgt-transforms.md",
            "Local-space operators" => "api/local-operators.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/lurlurlurrrrr/LurCGT.jl.git",
    devbranch = "main",
)
