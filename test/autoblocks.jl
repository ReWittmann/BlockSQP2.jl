using BlockSQP2
using BlockSQP2.NLPlayouts
using Optimization
import Optimization: SciMLBase
using SparseArrays
using Test

const ext = Base.get_extension(BlockSQP2, :OptimizationExtension)
const upper = sparse([1, 2], [3, 4], [0.0, 1.0], 4, 4)
const lower = sparse([3, 4], [1, 2], [0.0, 1.0], 4, 4)
const both = sparse([1, 2, 3, 4], [3, 4, 1, 2], [0.0, 1.0, 0.0, 1.0], 4, 4)
const detected = ([1, 3, 2, 4], Cint[0, 2, 4])
const fallback = (nothing, Cint[0, 4])

@testset "Stored Hessian structure" begin
    for pattern in (upper, lower, both, transpose(upper))
        @test ext.hessian_blocks(pattern, 4) == detected
    end
    isolated = sparse([1, 2], [4, 5], [0.0, 1.0], 6, 6)
    @test ext.hessian_blocks(isolated, 6) == ([1, 4, 2, 5, 3, 6], Cint[0, 2, 4, 5, 6])
    unsorted_bfs = sparse([1, 3], [4, 4], ones(2), 4, 4)
    @test ext.hessian_blocks(unsorted_bfs, 4) == ([1, 3, 4, 2], Cint[0, 3, 4])
    connected = sparse([2, 3, 4], [1, 1, 2], ones(3), 4, 4)
    @test ext.hessian_blocks(connected, 4) == (collect(1:4), Cint[0, 4])
    @test ext.hessian_blocks(spzeros(4, 4), 4) == (collect(1:4), Cint[0, 1, 2, 3, 4])
    @test_throws DimensionMismatch ext.hessian_blocks(spzeros(3, 4), 4)
end

const optimum = [1.0, 2.0, 3.0, 4.0]
const linear_term = [-2.0, 3.0, -7.0, 6.0]
const A = [1.0 1.0 0.0 0.0; 0.0 0.0 1.0 2.0]
const lb = [-2.0, 2.0, -4.0, -5.0]
const ub = [6.0, 7.0, 3.0, 9.0]
const x0 = [0.0, 3.0, 1.0, 5.0]

function quadratic(x, p)
    d = x - optimum
    return sum(abs2, d) + 0.5 * (d[1] * d[3] + d[2] * d[4]) + sum(linear_term .* d)
end

function gradient!(g, x, p)
    d = x - optimum
    g .= 2 .* d .+ 0.5 .* d[[3, 4, 1, 2]] .+ linear_term
    return nothing
end

function constrained_problem(; lag_pattern = both, obj_pattern = nothing)
    cons!(c, x, p) = (c .= A * x)
    jac!(J, x, p) = (J .= A)
    f = OptimizationFunction(quadratic; grad = gradient!, cons = cons!, cons_j = jac!,
        lag_hess_prototype = lag_pattern, hess_prototype = obj_pattern)
    return OptimizationProblem(f, x0; lb, ub, lcons = [3.0, 11.0], ucons = [3.0, 11.0])
end

@testset "Prototype and manual precedence" begin
    optimizer = BlockSQP2.Optimizer()
    problem = constrained_problem()
    structure(prob; kwargs...) = ext.automatic_structure(Optimization.init(prob, optimizer; kwargs...), 4, 2)
    @test structure(problem) == detected
    @test structure(constrained_problem(lag_pattern = nothing)) == fallback
    @test structure(constrained_problem(lag_pattern = nothing, obj_pattern = both)) == fallback
    @test structure(constrained_problem(lag_pattern = sparse(ones(4, 4)))) == fallback
    @test structure(problem; autoblocks = false) == fallback
    @test structure(problem; blockIdx = Cint[0, 1, 4]) == fallback
    @test structure(problem; vblocks = [BlockSQP2.vblock(4, false)]) == fallback
    hblock = BlockDescriptor{nlpHess}(tag = :variables)
    layout = NLPlayout((hblock,), to_Axis(TupleBD[(hblock, 4)]), (), to_Axis(TupleBD[]))
    @test structure(problem; layout) == fallback
    condenser = BlockSQP2.Condenser(
        [BlockSQP2.vblock(1, false), BlockSQP2.vblock(2, true), BlockSQP2.vblock(1, false)],
        [BlockSQP2.cblock(2)], [1, 3], [BlockSQP2.condensing_target(1, 0, 3, 0, 1)])
    @test structure(problem; condenser) == fallback
    @test_throws DimensionMismatch structure(constrained_problem(lag_pattern = spzeros(3, 4)))
end

@testset "Dense and sparse solver coordinates" begin
    problem = constrained_problem()
    optimizer = BlockSQP2.Optimizer()
    for sparse_mode in (false, true)
        baseline = solve(problem, optimizer; autoblocks = false,
            options = BlockSQP2.Options(sparse = sparse_mode, print_level = 0))
        iterates = Vector{Float64}[]
        function callback(state, loss)
            push!(iterates, copy(state.u))
            @test state.objective ≈ quadratic(state.u, nothing)
            @test loss ≈ state.objective
            return false
        end
        solution = solve(problem, optimizer; callback,
            options = BlockSQP2.Options(sparse = sparse_mode, print_level = 0))
        @test SciMLBase.successful_retcode(baseline)
        @test SciMLBase.successful_retcode(solution)
        @test solution.u ≈ optimum atol = 1e-6
        @test solution.u ≈ baseline.u atol = 1e-6
        @test solution.objective ≈ baseline.objective atol = 1e-6
        @test solution.objective ≈ 0.0 atol = 1e-6
        @test !isempty(iterates)
        @test last(iterates) ≈ optimum atol = 1e-6
        @test solution.original.bound_multiplier ≈ [0, -5, 10, 0] atol = 1e-5
        @test solution.original.multiplier ≈ [2, -3] atol = 1e-5
        @test solution.original.bound_multiplier ≈ baseline.original.bound_multiplier atol = 1e-5
        @test solution.original.multiplier ≈ baseline.original.multiplier atol = 1e-5

        empty!(iterates)
        stopped = solve(problem, optimizer;
            options = BlockSQP2.Options(sparse = sparse_mode, print_level = 0),
            callback = (state, loss) -> (push!(iterates, copy(state.u)); true))
        @test length(iterates) == 1
        @test stopped.u == only(iterates)
        @test stopped.objective ≈ quadratic(stopped.u, nothing)
    end
end

@testset "Unconstrained cache solve" begin
    f(x, p) = quadratic(x, p) - sum(linear_term .* (x - optimum))
    function grad!(g, x, p)
        gradient!(g, x, p)
        g .-= linear_term
    end
    problem = OptimizationProblem(OptimizationFunction(f; grad = grad!, hess_prototype = both), x0)
    cache = Optimization.init(problem, BlockSQP2.Optimizer())
    @test ext.automatic_structure(cache, 4, 0) == detected
    sol = Optimization.solve!(cache)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u ≈ optimum atol = 1e-6
end
