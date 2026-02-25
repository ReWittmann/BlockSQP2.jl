
using BlockSQP2
using Optimization
using ForwardDiff

_f(x,p) = x[1]^2 - 0.5*x[2]^2
_g(res, x,p) = res .=  [x[1] - x[2]]
x0 = Float64[10.0, 10.0]

optprob_wcons = OptimizationFunction(_f, AutoForwardDiff(), cons = _g)
prob = OptimizationProblem(optprob_wcons, x0, SciMLBase.NullParameters(), lcons = [0.0], ucons = [0.0])

options = BlockSQP2.Options(; sparse = true, max_conv_QPs = 1, indef_delay = 1, enable_linesearch = false)
sol_bsqp_wcons = solve(prob, BlockSQP2.Optimizer(); 
                       blockIdx=[0,1,2], 
                       options = options)
