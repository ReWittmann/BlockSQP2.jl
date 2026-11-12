module OptimizationExtension

using BlockSQP2, BlockSQP2.SparseArrays, BlockSQP2.NLPlayouts
using OptimizationBase, OptimizationBase.SciMLBase

SciMLBase.allowsbounds(::BlockSQP2Optimizer) = true
SciMLBase.allowsconstraints(::BlockSQP2Optimizer) = true
SciMLBase.allowscallback(::BlockSQP2Optimizer) = true
SciMLBase.requiresgradient(::BlockSQP2Optimizer) = true
SciMLBase.requiresconsjac(::BlockSQP2Optimizer) = true
SciMLBase.supports_opt_cache_interface(::BlockSQP2Optimizer) = true
SciMLBase.has_init(::BlockSQP2Optimizer) = true

# @info "Loading Optimization(Base).jl extension for BlockSQP2..."

function sparse_prototype(f, name)
    prototype = hasproperty(f, name) ? getproperty(f, name) : nothing
    return prototype isa AbstractMatrix && issparse(prototype) ? prototype : nothing
end

function hessian_blocks(pattern, n)
    size(pattern) == (n, n) ||
        throw(DimensionMismatch("Hessian prototype has size $(size(pattern)), expected ($n, $n)"))

    csc = pattern isa SparseArrays.AbstractSparseMatrixCSC ? pattern : sparse(pattern)
    neighbors = [Int[] for _ in 1:n]
    rows = rowvals(csc)
    for col in 1:n, k in nzrange(csc, col)
        row = rows[k]
        if row != col
            push!(neighbors[col], row)
            push!(neighbors[row], col)
        end
    end

    visited = falses(n)
    permutation = Int[]
    sizehint!(permutation, n)
    boundaries = Cint[0]
    for root in 1:n
        visited[root] && continue
        first_index = length(permutation) + 1
        push!(permutation, root)
        visited[root] = true
        head = first_index
        while head <= length(permutation)
            for neighbor in neighbors[permutation[head]]
                if !visited[neighbor]
                    visited[neighbor] = true
                    push!(permutation, neighbor)
                end
            end
            head += 1
        end
        sort!(view(permutation, first_index:length(permutation)))
        push!(boundaries, length(permutation))
    end
    return permutation, boundaries
end

function automatic_structure(cache, n, num_cons)
    args = cache.solver_args
    if args.autoblocks && isnothing(args.blockIdx) && isnothing(args.layout) &&
            isempty(args.vblocks) && isnothing(args.condenser)
        pattern = args._lag_hess_pattern
        if isnothing(pattern)
            pattern = sparse_prototype(cache.f, :lag_hess_prototype)
        end
        if isnothing(pattern) && num_cons == 0
            pattern = args._obj_hess_pattern
            if isnothing(pattern)
                pattern = sparse_prototype(cache.f, :hess_prototype)
            end
        end
        if !isnothing(pattern)
            permutation, boundaries = hessian_blocks(pattern, n)
            if length(boundaries) > 2
                identity_order = all(i -> permutation[i] == i, eachindex(permutation))
                return identity_order ? nothing : permutation, boundaries
            end
        end
    end
    return nothing, Cint[0, n]
end

function SciMLBase.__init(
            prob::SciMLBase.OptimizationProblem, opt::BlockSQP2Optimizer,
            ;
            callback = OptimizationBase.DEFAULT_CALLBACK,
            progress = false, maxiters = nothing,
            options::BlockSQP2.Options=BlockSQP2.Options(),
            blockIdx::Union{Vector{<:Integer}, Nothing} = nothing,
            vblocks::Union{Vector{BlockSQP2.vblock}, Nothing} = BlockSQP2.vblock[],
            condenser::Union{BlockSQP2.Condenser, Nothing} = nothing,
            layout::Union{NLPlayouts.NLPlayout, Nothing} = nothing,
            autoblocks::Bool = true,
            kwargs...
            )
    if isnothing(vblocks)
        vblocks = BlockSQP2.vblock[]
    end
    # Some OptimizationBase backends omit prototypes when instantiating prob.f.
    lag_pattern = sparse_prototype(prob.f, :lag_hess_prototype)
    obj_pattern = if isnothing(prob.ucons) || isempty(prob.ucons)
        sparse_prototype(prob.f, :hess_prototype)
    end
    return OptimizationCache(prob, opt;
        callback = callback, progress = progress, maxiters = maxiters,
        options=options, blockIdx = blockIdx, vblocks = vblocks, layout = layout,
        condenser = condenser, autoblocks = autoblocks,
        _lag_hess_pattern = lag_pattern, _obj_hess_pattern = obj_pattern, kwargs...)
end

function __map_optimizer_args!(cache::OptimizationCache,
    opt::BlockSQP2.Options;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    kwargs...)

    for j in kwargs
        if j.first in fieldnames(BlockSQP2.Options)
            setproperty!(opt, j.first, j.second)
        end
    end

    if !isnothing(maxiters)
        opt.maxiters = maxiters
    end

    if !isnothing(maxtime)
        @warn "common maxtime is currently not used by $(opt)"
    end

    if !isnothing(reltol)
        @warn "common reltol is currently not used by $(opt)"
    end

    if !isnothing(abstol)
        opt.opt_tol = abstol
    end

    return nothing
end


function SciMLBase.__solve(cache::OptimizationCache{<:BlockSQP2Optimizer})
    
    num_cons = cache.ucons === nothing ? 0 : length(cache.ucons)
    num_x = length(cache.u0)
    T = eltype(cache.u0)
    
    opts = cache.solver_args.options
    callback = cache.callback
    maxiters = OptimizationBase._check_and_convert_maxiters(cache.solver_args.maxiters)
    maxtime = OptimizationBase._check_and_convert_maxtime(cache.solver_args.maxtime)
    
    __map_optimizer_args!(cache, opts; callback = callback,
                maxiters = maxiters, maxtime = maxtime,
                abstol = cache.solver_args.abstol, reltol = cache.solver_args.reltol,
                cache.solver_args...)
    
    permutation, _blockIdx = automatic_structure(cache, num_x, num_cons)
    inverse_permutation = isnothing(permutation) ? nothing : invperm(permutation)
    to_solver(x) = isnothing(permutation) ? x : x[permutation]
    to_user(x) = isnothing(inverse_permutation) ? x : x[inverse_permutation]
    _vblocks = BlockSQP2.vblock[]
    _condenser = nothing
    
    if !isnothing(cache.solver_args.layout) 
        _layout = cache.solver_args.layout
        _blockIdx = hessBlockIndexZeroBased(_layout)
        _vblocks = create_vblocks(_layout)
        
        #Deactivate this for now, requiring explicit passing of a condenser.
        # _condenser = BlockSQP2.Condenser(_layout)
    end
    if !isnothing(cache.solver_args.condenser)
        _condenser = cache.solver_args.condenser
    end
    if length(cache.solver_args.vblocks) > 0
        _vblocks = cache.solver_args.vblocks
    end
    if !isnothing(cache.solver_args.blockIdx)
        _blockIdx = cache.solver_args.blockIdx
    end
    
    _loss = function (θ)
        x = cache.f(to_user(θ), cache.p)
        return only(x)
    end
    
    _g = function(θ)
        g_ = similar(θ)
        cache.f.grad(g_, to_user(θ))
        return to_solver(g_)
    end
    
    _cons = begin
        if num_cons > 0
            function(θ)
                c_ = zeros(eltype(θ), num_cons)
                cache.f.cons(c_, to_user(θ))
                return c_
            end
        else
            (x) -> zeros(eltype(x), 1) # dummy constraint that it doesnt crash, TODO: FIX THIS
        end
    end
    
    _jac_cons = begin
        if num_cons > 0
            function(θ)
                J = zeros(eltype(θ), num_cons, num_x)
                cache.f.cons_j(J, to_user(θ))
                return isnothing(permutation) ? J : J[:, permutation]
            end
        else
            (x) -> zeros(eltype(x), 1, num_x)
        end
    end
    
    u0 = to_solver(cache.u0)
    jac_g_row, jac_g_col, jac_row, jac_col = begin
        # if use_sparse_functions
        if opts.sparse
            #Hacky: Calculate constraint Jacobian for perturbed points to find structural nonzero elements
            u0_pert_1 = [x + 1e-6*rand() for x in u0]
            u0_pert_2 = [x + 1e-5*rand() for x in u0]
            u0_pert_3 = [x + 1e-4*rand() for x in u0]
            J_sparse = sparse(_jac_cons(u0) + _jac_cons(u0_pert_1) + _jac_cons(u0_pert_2) + _jac_cons(u0_pert_3))
            jac_row, jac_col, _ = findnz(J_sparse)
            J_sparse.rowval .- 1, J_sparse.colptr .- 1, jac_row, jac_col
        else
            Cint[], Cint[], nothing, nothing
        end
    end
    
    sparse_jac(x) = let jac_row = jac_row, jac_col = jac_col
        _J = _jac_cons(x)
        return [_J[i,j] for (i,j) in zip(jac_row,jac_col)]
    end
    
    num_cons = max(1, num_cons)
    _lb = isnothing(cache.lb) ? -Inf * ones(T, num_x) : cache.lb
    _ub = isnothing(cache.ub) ? Inf * ones(T, num_x) : cache.ub
    _lb_cons = isnothing(cache.lcons) ? -Inf * ones(T, num_cons) : cache.lcons
    _ub_cons = isnothing(cache.ucons) ? Inf * ones(T, num_cons) : cache.ucons
    
    _lambda_0 = zeros(T, num_x+num_cons)
    
    
    stats = BlockSQP2.Stats("./")
    
    _lb = to_solver(BlockSQP2.__lowerbounds(_lb))
    _ub = to_solver(BlockSQP2.__upperbounds(_ub))
    _u0 = BlockSQP2.__initial_values(u0)
    
    sqp_prob = BlockSQP2.Problem(_loss, _cons, _g, _jac_cons,
                            _lb, _ub, _lb_cons, _ub_cons,
                            _u0, _lambda_0; 
                            blockIdx = _blockIdx,
                            vblocks = _vblocks,
                            condenser = _condenser,
                            jac_g_row = jac_g_row, jac_g_colind = jac_g_col,
                            jac_g_nz = opts.sparse ? sparse_jac : BlockSQP2.fnothing
                            )
    
    meth = BlockSQP2.Solver(sqp_prob, opts, stats)
    
    BlockSQP2.init!(meth)
    t0 = time()
    if cache.callback == OptimizationBase.DEFAULT_CALLBACK
        ret = BlockSQP2.run!(meth, opts.maxiters, 1)
    else
        for i=1:opts.maxiters
            ret = BlockSQP2.run!(meth, 1, 1)
            _iterate = BlockSQP2.get_primal_solution(meth)
            _obj = _loss(_iterate)
            state = OptimizationBase.OptimizationState(iter = i,
                u = to_user(_iterate),
                objective = _obj)
            cret = cache.callback(state, _obj)
            (cret || BlockSQP2.is_success(ret)) && break
        end
    end
    BlockSQP2.finish!(meth)
    t1 = time()
    
    x_opt = BlockSQP2.get_primal_solution(meth)
    lambda = BlockSQP2.get_dual_solution_full(meth)
    bound_multiplier = to_user(lambda[1:num_x])
    multiplier = lambda[(num_x + 1):end]
    f_opt = _loss(x_opt)
    retcode = BlockSQP2.is_success(ret) ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.Default
    
    original = (ret = ret, multiplier = multiplier, solve_time = t1 - t0,
        solve_it = Int64(BlockSQP2.get_itCount(meth)), bound_multiplier = bound_multiplier)
    return SciMLBase.build_solution(cache, cache.opt, to_user(x_opt), f_opt; original, retcode)
end

end # module OptimizationExtension
