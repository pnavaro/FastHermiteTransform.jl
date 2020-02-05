using Documenter
using FastHermiteTransform

makedocs(
    sitename = "FastHermiteTransform",
    format = Documenter.HTML(),
    modules = [FastHermiteTransform]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
