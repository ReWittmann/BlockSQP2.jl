**BlockSQP2.jl** -- A Julia interface to [blockSQP2](https://github.com/ReWittmann/blockSQP2), a nonlinear  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; programming solver based on [blockSQP](https://github.com/djanka2/blockSQP).  
Copyright (c) 2025 Reinhold Wittmann <reinhold.wittmann@ovgu.de> and Christoph Plate <christoph.plate@ovgu.de>

## Installation

This package requires the compiled binary of the C interface to [blockSQP2](https://github.com/ReWittmann/blockSQP2). One can either  
 - Fetch this package as a submodule of [blockSQP2](https://github.com/ReWittmann/blockSQP2), letting the build-process place the binary in BlockSQP2.jl/bin  
 - Depend on [blockSQP2_jll](https://github.com/JuliaRegistries/General/tree/master/jll/B/blockSQP2_jll)
  
The binary may link to, or include, compiled code subject to distinct licenses, see blockSQP2/README.md for a list of dependencies and their licenses.

Afterwards, the package can be managed through the Julia package manager as usual, see <https://docs.julialang.org/en/v1/stdlib/Pkg/>.

### Licensing
BlockSQP2.jl is published under the very permissive zlib license, see LICENSE.txt.

### Example from blockSQP
**min** &nbsp; x<sub>1</sub><sup>2</sup> - 0.5x<sub>2</sub><sup>2</sup>  
&nbsp; **s.t.** &nbsp; -inf < x<sub>1</sub>, x<sub>2</sub> < inf  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 0 ≤ x<sub>1</sub> - x<sub>2</sub> ≤ 0
```julia
using BlockSQP2

_f(x,p) = x[1]^2 - 0.5*x[2]^2
_g(res,x,p) = res .=  [x[1] - x[2]]
x0 = Float64[10.0, 10.0]
```

BlockSQP2.jl can be called either via the [Optimization.jl](https://docs.sciml.ai/Optimization/stable/) interface ...

```julia
using Optimization
using ForwardDiff

optprob = OptimizationFunction(_f, AutoForwardDiff(), cons = _g)
prob = OptimizationProblem(optprob, x0, SciMLBase.NullParameters(), lcons = [0.0], ucons = [0.0])

options = BlockSQP2.Options(opt_tol = 1.0e-12, 
                            feas_tol = 1.0e-12, 
                            max_conv_QPs = 1, 
                            indef_delay = 1, 
                            enable_linesearch = false)

sol_bsqp = solve(prob, BlockSQP2.Optimizer(); 
                 blockIdx = [0,1,2], # block indices of the block-diagonal Hessian
                 options = options)
```
... or directly
```julia
f(x) = _f(x, nothing)
g(x) = (res = [0.]; _g(res, x, nothing); res)

# Derivatives
grad_f = x::Array{Float64, 1} -> Float64[2*x[1], -x[2]]
jac_g = x::Array{Float64, 1} -> Float64[1 -1]

## Sparse constraint Jacobian in CCS format
jac_g_nz = x::Array{Float64, 1} -> Float64[1, -1]
jac_g_row = Int32[0, 0]
jac_g_colind = Int32[0, 1, 2]

# Bounds
lb_var = Float64[-Inf, -Inf]
ub_var = Float64[Inf, Inf]
lb_con = Float64[0.0]
ub_con = Float64[0.0]

lambda0 = Float64[0., 0., 0.]

prob = BlockSQP2.Problem(f, g, grad_f, jac_g, 
                         lb_var, ub_var, lb_con, ub_con, 
                         x0, lambda0; 
                         jac_g_nz = jac_g_nz, 
                         jac_g_row = jac_g_row, 
                         jac_g_colind = jac_g_colind, 
                         blockIdx = Int32[0, 1, 2])

options = BlockSQP2.Options(opt_tol = 1.0e-12, 
                            feas_tol = 1.0e-12, 
                            max_conv_QPs = 1, 
                            indef_delay = 1, 
                            enable_linesearch = false)


stats = BlockSQP2.Stats("./")

meth = BlockSQP2.Solver(prob, options, stats)

init!(meth)
ret = run!(meth, Int32(100), Int32(1))
finish!(meth)

x_opt = get_primal_solution(meth)
lam_opt = get_dual_solution(meth)

```
### See also
[Corleone.jl](https://github.com/SciML/Corleone.jl) - setup optimal control problems, apply multiple-shooting to obtain structured nonlinear optimization problems, and much more.  
  
*/examples/lotka_oc_Corleone.jl* - example of using Corleone 0.0.3 to setup and parameterize an optimal control problem via [multiple shooting](https://www.sciencedirect.com/science/article/pii/S1474667017612059), and solving it efficiently with BlockSQP2.  
  
[blockSQP2 paper preprint](https://optimization-online.org/?p=32379)
