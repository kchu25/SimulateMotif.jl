using SimulateMotif
using Documenter

DocMeta.setdocmeta!(SimulateMotif, :DocTestSetup, :(using SimulateMotif); recursive=true)

makedocs(;
    modules=[SimulateMotif],
    authors="Shane Kuei-Hsien Chu",
    repo="https://github.com/kchu25/SimulateMotif.jl/blob/{commit}{path}#{line}",
    sitename="SimulateMotif.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/SimulateMotif.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/SimulateMotif.jl",
    devbranch="main",
)
