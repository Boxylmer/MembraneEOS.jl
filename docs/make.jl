using Documenter
using MembraneEOS

makedocs(;
    modules=[MembraneBase],
    authors="Will <william.joseph.box@gmail.com> and contributors",
    repo="https://github.com/Boxylmer/MembraneEOS.jl/blob/{commit}{path}#{line}",
    sitename="MembraneBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Boxylmer.github.io/MembraneEOS.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Boxylmer/MembraneEOS.jl.git",
    devbranch="master",
    
)