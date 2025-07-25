using Pkg;
# Activate current blockSQP environment
Pkg.activate(@__DIR__)
# Add binary
Pkg.develop(url="https://github.com/chplate/qpoases_mumps.git")
Pkg.develop(url="https://github.com/chplate/blocksqp_mumps.git")
