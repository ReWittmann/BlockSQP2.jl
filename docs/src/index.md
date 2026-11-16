# BlockSQP2.jl

BlockSQP2.jl provides Julia callbacks and an Optimization.jl interface to the
native blockSQP2 nonlinear programming solver. It solves problems of the form

```math
\min_x f(x) \quad\text{subject to}\quad
\ell_x \le x \le u_x,\qquad \ell_g \le g(x) \le u_g.
```

The solver uses sequential quadratic programming (SQP), with blockwise Hessian
approximations. For multiple shooting problems, structural metadata enables
condensing of the QP subproblems.

## Installation

This development checkout requires Julia 1.11 or later and a locally built native
interface in `BlockSQP2.jl/bin`. Follow `docs/README.md` for native and documentation
build commands. From the repository root, add the local package to an environment:

```julia
using Pkg
Pkg.develop(path="BlockSQP2.jl")
```

The source currently selects the development library loader. Its commented
release loader uses `blockSQP2_jll`; switching loaders is a source configuration
change.

## A direct solve

This two-variable example minimizes ``x_1^2 - x_2^2/2`` with ``x_1=x_2``.
The direct interface accepts allocating Julia callbacks. Dense Jacobians need
`sparse=false`, since the default options select sparse evaluation.

```@example getting_started
using BlockSQP2

f(x) = x[1]^2 - 0.5*x[2]^2
g(x) = [x[1] - x[2]]
grad_f(x) = [2*x[1], -x[2]]
jac_g(x) = reshape([1.0, -1.0], 1, 2)

prob = BlockSQP2.Problem(
    f, g, grad_f, jac_g,
    fill(-Inf, 2), fill(Inf, 2), [0.0], [0.0],
    [10.0, 10.0], zeros(3); blockIdx=[0, 1, 2],
)
opts = BlockSQP2.Options(sparse=false, print_level=0, print_colored=false)
x_opt = mktempdir() do directory
    stats = BlockSQP2.Stats(joinpath(directory, ""))
    solver = BlockSQP2.Solver(prob, opts, stats)
    init!(solver)
    result = run!(solver, 100, 0)
    finish!(solver)
    @assert BlockSQP2.is_success(result)
    get_primal_solution(solver)
end
@assert maximum(abs, x_opt) < 1e-4
round.(x_opt; digits=4)
```

Construct settings before creating the solver: `Solver` transfers them to native
objects. `init!` initializes the iteration, `run!` returns a status, and `finish!`
finalizes reporting. Julia finalizers release native allocations. Solution
accessors return copies; the full dual accessor puts variable-bound multipliers
before general-constraint multipliers and reverses the native sign.

## Optimization.jl

Loading Optimization.jl activates an extension for [`BlockSQP2.Optimizer`](@ref).
For an existing `OptimizationProblem`, use:

```julia
using Optimization, ForwardDiff, BlockSQP2
fun = OptimizationFunction((x, p) -> sum(abs2, x), AutoForwardDiff())
problem = OptimizationProblem(fun, [1.0, 2.0])
solution = solve(problem, BlockSQP2.Optimizer(); maxiters=100)
```

The extension accepts `options`, `blockIdx`, `vblocks`, `condenser`, and `layout`
keywords, supports callbacks and the cache interface, and obtains derivatives
through Optimization.jl. `abstol` maps to the optimality tolerance; `reltol` and
`maxtime` are currently unused. See the package tests for constrained examples.

Continue with the [Lotka–Volterra walkthrough](@ref), or consult the
[API reference](@ref).
