using ScalarWave1D
using Documenter

DocMeta.setdocmeta!(ScalarWave1D, :DocTestSetup, :(using ScalarWave1D); recursive=true)

makedocs(;
    modules=[ScalarWave1D],
    authors="Stamatis Vretinaris",
    repo="https://github.com/svretina/ScalarWave1D.jl/blob/{commit}{path}#{line}",
    sitename="ScalarWave1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://svretina.github.io/ScalarWave1D.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/svretina/ScalarWave1D.jl",
    devbranch="master",
)
