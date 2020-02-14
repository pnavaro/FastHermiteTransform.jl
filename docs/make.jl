using Documenter
using FastHermiteTransform

makedocs(
    sitename = "FastHermiteTransform.jl",
    format = Documenter.HTML(),
    modules = [FastHermiteTransform],
    pages = ["Documentation"    => "index.md"]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pnavaro/FastHermiteTransform.jl.git",
)
