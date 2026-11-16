# Lotka–Volterra walkthrough

The supplied `examples/example_lotka.jl` solves an optimal harvesting problem by
multiple shooting. This page follows its model, layout, derivative setup, and
solver calls. The runnable `docs/lotka.jl` contains the same construction with
the plotting and comparison code removed. Run it with the documentation
environment as described in `docs/README.md`.

## Model and integration

Two populations start at ``x(0)=(0.5,0.7)``. A scalar harvesting control
``0\le u\le1`` acts over ``0\le t\le12``:

```math
\dot x_1=x_1-x_1x_2-0.4ux_1,\qquad
\dot x_2=-x_2+x_1x_2-0.2ux_2.
```

The objective penalizes distance from the equilibrium population ``(1,1)``:

```math
J=\int_0^{12} ((x_1(t)-1)^2+(x_2(t)-1)^2)\,dt.
```

`lotka_rhs(x,u)` evaluates the dynamics, and `lotka_quad(x,u)` returns a
one-element vector containing the running cost. `RK4_step` advances both with
classical fourth-order Runge–Kutta; `RK4_M` takes two substeps per shooting
interval. `ODEsol_multi` integrates the intervals independently, using the
candidate state at each interval's left endpoint. Its threaded loop can use
multiple Julia threads. The control is constant within each interval.

## Multiple shooting variables and constraints

With `N=100`, there are 100 controls and 100 two-component state vectors. The
initial state is fixed data, so the decision vector has 300 entries:

```math
z=(u_0,\;x_1,u_1,\;x_2,u_2,\;\ldots,\;x_{99},u_{99},\;x_{100}).
```

Here a subscript denotes a shooting node, while each ``x_k`` contains both
populations. Define ``\Phi(x_k,u_k)`` as the RK4 interval map. The 200 equality
constraints enforce

```math
x_{k+1}-\Phi(x_k,u_k)=0,\qquad k=0,\ldots,99.
```

`objective` reshapes the flat state and control slices, prepends the fixed
initial state, integrates all intervals, and sums their quadrature values.
`shooting_constraints` flattens the state defects in Julia's column order.
All decision variables have lower bound zero; controls have upper bound one
and states have upper bound infinity. The starting guess repeats the initial
population at every state node and uses zero controls. This guess need not
satisfy the matching constraints.

## Describe the structure with NLPlayouts

The layout describes the existing callbacks' vector ordering. It does not
generate the objective or constraints. `BlockDescriptor`s identify nested
sections and carry links used for structural checks and condensing:

```julia
using BlockSQP2, BlockSQP2.NLPlayouts
matchings = BlockDescriptor{nlpMatchings}(tag=:matchings)
MSsys = BlockDescriptor{nlpMultipleShootingDF}(matchings=matchings, tag=:MSsys)
h0 = BlockDescriptor{nlpHess}(parent=MSsys, tag=:h₀)
u0 = BlockDescriptor{nlpMSfree}(parent=h0, tag=:u₀)
c1 = BlockDescriptor{nlpMatching}(parent=matchings, input=[u0], tag=:m₁)
```

Each Hessian block groups the variables used by one interval cost and map:
first ``u_0``, then ``(x_k,u_k)``, and finally ``x_{100}``. The terminal block
is retained even though this model has no terminal cost. Controls are free
blocks; states are dependent blocks, linked to the preceding matching block
with `matching=ck`. A matching block's `input` lists the variable descriptors
on which its interval map depends. At the first interval the initial state is
fixed, so only `u0` is listed.

The example assembles nested `(descriptor, size_or_children)` pairs in
`TupleBD[]` pre-layouts. The tree and each descriptor's `parent` must agree.
After the loop constructs all stages:

```julia
vLayout = to_Axis(vPreLayout)
cLayout = to_Axis(cPreLayout)
layout = NLPlayout(
    (blockDescriptors(vPreLayout)...,), vLayout,
    (blockDescriptors(cPreLayout)...,), cLayout,
)
```

`axsubrange(vLayout, descriptor)` returns ordinary **one based Julia indices**
into the flat vector. `ComponentArray(data, vLayout)` attaches names to that
vector, enabling `x_opt[state]` and `view(x_start, state)`.
`hessBlockIndexZeroBased(layout)` produces the boundary array for the native
solver: `[0, 1, 4, …, 298, 300]`, describing 101 Hessian blocks.

## Differentiate the callbacks

The objective gradient uses ForwardDiff. The constraint Jacobian uses sparse
forward differentiation with traced sparsity and column coloring:

```julia
using ForwardDiff, SparseConnectivityTracer, SparseMatrixColorings
using DifferentiationInterface, SparseArrays

grad_f = x -> ForwardDiff.gradient(objective, x)
sparse_forward_backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector=TracerSparsityDetector(),
    coloring_algorithm=GreedyColoringAlgorithm(),
)
jac_g0 = jacobian(shooting_constraints, sparse_forward_backend, x_start)
ROW = jac_g0.rowval .- 1
COLIND = jac_g0.colptr .- 1
jac_gNZ(x) = jacobian(shooting_constraints, sparse_forward_backend, x).nzval
```

The native interface uses compressed sparse column (CSC) storage. Subtract one
from Julia's row indices and column pointers, but keep the values in the same
order. The sparsity pattern must stay fixed for the whole solve, including
entries whose values happen to vanish. `jac_gNZ` must return exactly the number
of entries described by `ROW`. No exact Hessian callback is needed for the
default SR1 approximation and BFGS fallback.

## Condense and solve

`create_vblocks(layout)` converts leaf descriptors to the native variable block
metadata. Full condensing eliminates dependent state steps within each target;
partial condensing keeps a shorter shooting structure. The example uses:

```julia
condenser = BlockSQP2.PartialCondenser(layout, 4) # four regions per target
prob = BlockSQP2.Problem(
    objective, shooting_constraints, grad_f, BlockSQP2.fnothing,
    collect(lb_var), collect(ub_var), lb_con, ub_con,
    collect(x_start), zeros(axlength(vLayout) + axlength(cLayout));
    blockIdx=hessBlockIndexZeroBased(layout),
    jac_g_row=ROW, jac_g_colind=COLIND, jac_g_nz=jac_gNZ,
    vblocks=BlockSQP2.create_vblocks(layout), condenser=condenser,
)
opts = BlockSQP2.sparse_options()
opts.max_conv_QPs = 4
opts.conv_strategy = "reduced_regularization"
```

The dense Jacobian argument is unused on this sparse path, hence the placeholder
`BlockSQP2.fnothing`. Condensers built from layouts currently emit an experimental
feature warning. Set `condenser=nothing` to solve the original sparse QPs.
The original example additionally enables automatic scaling and builds
zero based variable groups with `primitive_get_vsystem`. The headless script
uses the default `automatic_scaling=false` and omits that optional grouping.

```julia
stats = BlockSQP2.Stats("./")
meth = BlockSQP2.Solver(prob, opts, stats)
init!(meth)
result = run!(meth, 200, 0)
finish!(meth)
@assert BlockSQP2.is_success(result)
x_opt = ComponentArray(get_primal_solution(meth), layout.vLayout)
```

The state trajectory includes the fixed initial value; the control trajectory
has one value per interval:

```julia
opt_x = hcat(lotka_vparams[:x_init], (x_opt[st] for st in states)...)
opt_u = hcat((x_opt[ctrl] for ctrl in controls)...)
Tgrid = range(0, 12; length=N+1)
```

Check the returned status, the matching residuals, and the objective before
interpreting the trajectory. The existing Lotka test checks an objective near
`1.344408` with absolute tolerance `1e-5` for this discretization. The headless
script also checks bounds and matching feasibility.

## Plotting and comparing solvers

The original example solves the same NLP with Ipopt through Optimization.jl
before running BlockSQP2. Its Ipopt path uses dense Jacobians, so the printed
runtimes do not provide an equivalent sparse comparison.

At the end, the example loads CairoMakie and GLMakie and waits for a display
window. For a headless plot, load only CairoMakie, plot the populations with
`lines!`, plot interval controls with `stairs!`, and call `save("lotka.png", fig)`.
The documentation environment intentionally contains only the packages needed
for building the site and running the headless solve; plotting and the Ipopt
comparison use `examples/Project.toml`.
