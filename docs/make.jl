using Documenter
using DiscreteChoice

makedocs(
    format = :html,
    sitename = "DiscreteChoice.jl",
    pages = [
        "index.md",
        "getting_started.md"
    ]
)

deploydocs(
    repo = "github.com/brilhana/DiscreteChoice.jl.git",
    julia  = "1.0",
    latest = "master",
    target = "build",
    deps = nothing,
    make = nothing
)