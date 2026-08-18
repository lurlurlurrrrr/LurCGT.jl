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
        "Public API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/lurlurlurrrrr/LurCGT.jl.git",
    devbranch = "main",
)
