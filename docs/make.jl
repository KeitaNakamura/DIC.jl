using Documenter
using DIC

# Setup for doctests in docstrings
DocMeta.setdocmeta!(DIC, :DocTestSetup, recursive = true,
    quote
        using DIC
    end
)

makedocs(;
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [DIC],
    sitename = "DIC.jl",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "DIC Searching.md"
            "Utilities.md"
        ],
    ],
    doctest = true, # :fix
)

deploydocs(
    repo = "github.com/KeitaNakamura/DIC.jl.git",
    devbranch = "main",
)
