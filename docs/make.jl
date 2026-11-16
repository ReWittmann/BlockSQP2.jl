using BlockSQP2
using Documenter

makedocs(;
    root = @__DIR__,
    sitename = "BlockSQP2.jl",
    modules = [BlockSQP2, BlockSQP2.NLPlayouts],
    remotes = nothing,
    checkdocs = :all,
    format = Documenter.HTML(
        prettyurls = false,
        edit_link = nothing,
        repolink = "https://github.com/ReWittmann/blockSQP2",
    ),
    pages = [
        "Getting started" => "index.md",
        "Lotka–Volterra walkthrough" => "lotka.md",
        "API reference" => "api.md",
    ],
)
