# Headless adaptation of examples/example_lotka.jl. See docs/src/lotka.md.
using BlockSQP2
using BlockSQP2.NLPlayouts

const lotka_params = Dict{Symbol, Float64}(
    :c0 => 0.4,
    :c1 => 0.2,
    :t0 => 0.,
    :t1 => 12.
)

const lotka_vparams = Dict{Symbol, Vector{Float64}}(
    :x_init => [0.5, 0.7]
)

function lotka_rhs(x,u)
    c0, c1 = (lotka_params[x] for x in (:c0, :c1))
    return [x[1] - x[1]*x[2] - u[1]*c0*x[1], -x[2] + x[1]*x[2] - u[1]*c1*x[2]]
end
const lotka_nx = 2
const lotka_nu = 1

function lotka_quad(x,u)
    return [(x[1] - 1)^2 + (x[2] - 1)^2]
end

function RK4_step(rhs, quad, x0, u, DT = 1.0)
    K1, K1q = rhs(x0, u), quad(x0, u)
    K2, K2q = rhs(x0 + DT/2 * K1, u), quad(x0 + DT/2 * K1, u)
    K3, K3q = rhs(x0 + DT/2 * K2, u), quad(x0 + DT/2 * K2, u)
    K4, K4q = rhs(x0 + DT * K3, u), quad(x0 + DT * K3, u)
    x1 = x0 + DT/6*(K1 + 2*K2 + 2*K3 + K4)
    Q  = DT/6*(K1q + 2*K2q + 2*K3q + K4q)
    return x1, Q
end

function RK4_M(rhs, quad, x0, u, DT = 1.0, M = 2)
    x, q = RK4_step(rhs, quad, x0, u, DT/M)
    for i = 2:M
        x, _q = RK4_step(rhs, quad, x, u, DT/M)
        q += _q
    end
    return x, q
end

function ODEsol_multi(rhs, quad, xk::AbstractArray{DTP}, uk::AbstractArray{DTP}, T = 12.0, M = 2) where DTP <: Number
    _N_stages = size(xk, 2)
    DT = T/_N_stages
    TMP = Vector{Tuple{Vector{DTP}, Vector{DTP}}}(undef, _N_stages)
    Threads.@threads for i in 1:_N_stages
        TMP[i] = RK4_M(rhs, quad, xk[:,i], uk[:,i], DT, M)
    end
    
    q_out = sum(last.(TMP))
    x_out = reduce(hcat, first.(TMP))
    return x_out, q_out
end

const N = 100
const nx = 2
const nu = 1
const nq = 1
_subscript(i::Integer) = (i |> digits |> reverse .|> dgt->Char(0x2080+dgt)) |> join

vPreLayout = TupleBD[]
cPreLayout = TupleBD[]

states = BlockDescriptor[]
controls = BlockDescriptor[]

matchings = BlockDescriptor{nlpMatchings}(tag = :matchings)
MSsys = BlockDescriptor{nlpMultipleShootingDF}(matchings = matchings, tag = :MSsys)

h0 = BlockDescriptor{nlpHess}(parent = MSsys, tag = :h₀)
u0 = BlockDescriptor{nlpMSfree}(parent = h0, tag = :u₀)

push!(vPreLayout, (h0, [(u0, nu)]))
push!(controls, u0)

c1 = BlockDescriptor{nlpMatching}(tag = :m₁, parent = matchings, input = [u0])
push!(cPreLayout, (c1, nx))

hk, uk, ck = h0, u0, c1
xk = BlockDescriptor()
for i = 1:(N-1)
    global hk = BlockDescriptor{nlpHess}(parent = MSsys, tag = Symbol(:h, _subscript(i)))
    global xk = BlockDescriptor{nlpMSdependent}(parent = hk, matching = ck, tag = Symbol(:x, _subscript(i)))
    global uk = BlockDescriptor{nlpMSfree}(parent = hk, tag = Symbol(:u, _subscript(i)))
    
    push!(vPreLayout, (hk, [(xk, nx), (uk, nu)]))
    
    push!(states, xk)
    push!(controls, uk)
    
    global ck = BlockDescriptor{nlpMatching}(parent = matchings, input = [xk, uk], tag = Symbol(:m, _subscript(i+1)))
    push!(cPreLayout, (ck, nx))
end

hN = BlockDescriptor{nlpHess}(parent = MSsys, tag = Symbol(:h, _subscript(N)))
xN = BlockDescriptor{nlpMSdependent}(parent = hN, matching = ck, tag = Symbol(:x, _subscript(N)))
# uN = BlockDescriptor{nlpMSfree}(parent = hN, tag = Symbol(:u, _subscript(N)))
push!(vPreLayout, (hN, [(xN, nx)]))#, (uN, 0)]))
push!(states, xN)

vPreLayout = [(MSsys, vPreLayout)]
cPreLayout = [(matchings, cPreLayout)]

vLayout = to_Axis(vPreLayout)
cLayout = to_Axis(cPreLayout)
vBlocks = blockDescriptors(vPreLayout)
cBlocks = blockDescriptors(cPreLayout)

layout = NLPlayout((vBlocks...,), vLayout, (cBlocks...,), cLayout)
nVar = axlength(vLayout)
lotka_MS(x,u) = ODEsol_multi(lotka_rhs, lotka_quad, x, u, lotka_params[:t1] - lotka_params[:t0], 2)
stateidx = splat(vcat)(states .|> st->collect(axsubrange(vLayout, st)))
controlidx = splat(vcat)(controls .|> ctrl->collect(axsubrange(vLayout, ctrl)))


function objective(arg_x)
    let _sidx = stateidx, _cidx = controlidx
        x0 = lotka_vparams[:x_init]
        x = hcat(x0, reshape(arg_x[_sidx], 2, :))
        u = reshape(arg_x[_cidx], 1, :)
        
        _, Q = lotka_MS(x[:,1:end-1], u)
        return first(Q)
    end
end

function shooting_constraints(arg_x)
    let _sidx = stateidx, _cidx = controlidx
        x0 = lotka_vparams[:x_init]
        x = hcat(x0, reshape(arg_x[_sidx], 2, :))
        u = reshape(arg_x[_cidx], 1, :)
        
        S, _ = lotka_MS(x[:, 1:end-1], u)
        return collect(Iterators.flatten(x[:,2:end] - S))
    end
end

x_init = lotka_vparams[:x_init]
x_start = ComponentArray(zeros(axlength(vLayout)), vLayout)
# u_start = zeros(1,N)

lb_var = ComponentArray(repeat([0.], axlength(vLayout)), vLayout)
ub_var = ComponentArray(repeat([Inf], axlength(vLayout)), vLayout)
for state in states
    view(x_start, state)[:] = x_init
end
for control in controls
    view(ub_var, control)[:] .= 1
end

lb_con, ub_con = zeros(axlength(cLayout)), zeros(axlength(cLayout))


using ForwardDiff

f = objective
g = shooting_constraints
grad_f = x->ForwardDiff.gradient(f, x)


using SparseConnectivityTracer
using SparseMatrixColorings
using DifferentiationInterface
using SparseArrays

sparse_forward_backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm()
)


jac_g0 = jacobian(g, sparse_forward_backend, x_start)
ROW = jac_g0.rowval .-1
COLIND = jac_g0.colptr .-1

jac_gNZ(x) = jacobian(g, sparse_forward_backend, x).nzval


# Condense the 100 shooting intervals into four target regions.
condenser = BlockSQP2.PartialCondenser(layout, 4)
prob = BlockSQP2.Problem(
    f, g, grad_f, BlockSQP2.fnothing,
    collect(lb_var), collect(ub_var), lb_con, ub_con,
    collect(x_start), zeros(axlength(vLayout) + axlength(cLayout));
    blockIdx=hessBlockIndexZeroBased(layout),
    jac_g_row=ROW, jac_g_colind=COLIND, jac_g_nz=jac_gNZ,
    vblocks=BlockSQP2.create_vblocks(layout), condenser=condenser,
)
opts = BlockSQP2.sparse_options()
opts.max_conv_QPs = 4
opts.conv_strategy = "reduced_regularization"
opts.print_level = 0
opts.print_colored = false

# Keep native solver output out of the source tree.
x_opt = mktempdir() do directory
    stats = BlockSQP2.Stats(joinpath(directory, ""))
    meth = BlockSQP2.Solver(prob, opts, stats)
    init!(meth)
    result = run!(meth, 200, 0)
    finish!(meth)
    @assert BlockSQP2.is_success(result) result
    get_primal_solution(meth)
end
@assert isapprox(f(x_opt), 1.344408; atol=1e-5)
@assert maximum(abs, g(x_opt)) < 1e-5
@assert all(collect(lb_var) .- 1e-5 .<= x_opt)
@assert all(x_opt .<= collect(ub_var) .+ 1e-5)
x_opt = ComponentArray(x_opt, layout.vLayout)
opt_x = hcat(lotka_vparams[:x_init], (x_opt[st] for st in states)...)
opt_u = hcat((x_opt[ctrl] for ctrl in controls)...)
println("Lotka objective: ", f(x_opt))
println("Maximum matching residual: ", maximum(abs, g(x_opt)))
