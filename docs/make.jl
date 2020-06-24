using NetMSA
using Documenter

makedocs(;
    modules=[NetMSA],
    authors="Aadam <aadimator@gmail.com> and contributors",
    repo="https://github.com/aadimator/NetMSA.jl/blob/{commit}{path}#L{line}",
    sitename="NetMSA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aadimator.github.io/NetMSA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aadimator/NetMSA.jl",
)
