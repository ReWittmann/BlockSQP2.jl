using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(path = joinpath(@__DIR__, ".."))

Pkg.add(name = "Corleone", version="0.0.3")
Pkg.add("LuxCore")
Pkg.add("OrdinaryDiffEq")
Pkg.add("SciMLSensitivity")
Pkg.add("ComponentArrays")


Pkg.add("Optimization")
Pkg.add("OptimizationMOI")
Pkg.add("Ipopt")

Pkg.add("QPALM")

Pkg.add("DifferentiationInterface")
Pkg.add("SparseConnectivityTracer")
Pkg.add("SparseMatrixColorings")
Pkg.add("SparseArrays")
Pkg.add("ForwardDiff")
Pkg.add("Zygote")

Pkg.add("CairoMakie")
