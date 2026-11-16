"""
    Problem(f, g, grad_f, jac_g, lb_var, ub_var, lb_con, ub_con, x0, lambda0; kwargs...)

Describe `min f(x)` subject to variable bounds and `lb_con ≤ g(x) ≤ ub_con`.
Callbacks return a scalar objective, constraint vector, gradient vector, and an
`nCon × nVar` Jacobian. Use `Float64` vectors for bounds and initial values;
`lambda0` has `nVar + nCon` entries (variable bounds first).

`blockIdx` gives contiguous Hessian block boundaries, normally zero based with
endpoints `0` and `nVar`; one based boundaries are converted automatically.
For sparse evaluation, supply `jac_g_nz(x)` and fixed, zero based CSC arrays
`jac_g_row` and `jac_g_colind`, and use `Options(sparse=true)`.
Otherwise supply `jac_g(x)` and use `Options(sparse=false)`.
Optional `vblocks`, `condenser`, and zero based `vgroup_ind` describe structure
for condensing and scaling. Also available as `BlockSQP2Problem`.
"""
mutable struct Problem
    
    nVar::Cint
    nCon::Cint
    nnz::Cint
    blockIdx::Vector{Cint}
    
    vblocks::Vector{vblock}
    condenser::Union{Condenser, Nothing}
    
    vgroup_types::Vector{String}
    vgroup_ind::Vector{Vector{Cint}}
    
    lb_var::Vector{Cdouble}
    ub_var::Vector{Cdouble}
    lb_con::Vector{Cdouble}
    ub_con::Vector{Cdouble}
    lb_obj::Cdouble
    ub_obj::Cdouble
    
    f::Function #Vector{Cdouble}[nVar] -> Cdouble
    g::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nCon]
    grad_f::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nVar]
    jac_g::Function #Vector{Cdouble}[nVar] -> Matrix{Cdouble, 2}[nCon, nVar]
    last_hessBlock::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[lastBlocksize*(lastBlocksize + 1)/2] #lower diagonal elements
    hess::Function #Vector{Cdouble}[nVar] -> Vector{Vector{Cdouble}}[blocksize*(blocksize + 1)/2][length(blockIdx) - 1] #See utils.jl lower_to_full!, full_to_lower!
    
    # constrVioReducer::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nVar]
    
    jac_g_nz::Function
    jac_g_row::Vector{Cint}
    jac_g_colind::Vector{Cint}
    
    x0::Vector{Cdouble}
    lambda0::Vector{Cdouble}
    
    function Problem(f::Function,
                     g::Function,
                     grad_f::Function,
                     jac_g::Function,
                     lb_var::Vector{FLOAT_T},
                     ub_var::Vector{FLOAT_T},
                     lb_con::Vector{FLOAT_T},
                     ub_con::Vector{FLOAT_T},
                     x0::Vector{FLOAT_T},
                     lambda0::Vector{FLOAT_T};
                     lb_obj::AbstractFloat = -Inf, 
                     ub_obj::AbstractFloat = Inf,
                     blockIdx::Vector{INT_T_1} = Int64[0, length(lb_var)],
                     vblocks::Vector{vblock} = vblock[],
                     condenser::Union{Condenser, Nothing} = nothing,
                     vgroup_types::Vector{String} = String[],
                     vgroup_ind::Vector{Vector{INT_T_2}} = Vector{Int64}[],
                     jac_g_row::Vector{INT_T_3} = Int64[],
                     jac_g_colind::Vector{INT_T_4} = Int64[],
                     jac_g_nz::Function = fnothing, 
                    #  constrVioReducer::Function = fnothing,
                     last_hessBlock::Function = fnothing, 
                     hess::Function = fnothing
                     ) where {FLOAT_T <: AbstractFloat, INT_T_1 <: Integer, INT_T_2 <: Integer, INT_T_3 <: Integer, INT_T_4 <: Integer}
                        if blockIdx[1] == 1
                            blockIdx = copy(blockIdx) .-1
                        end
                        @assert blockIdx[1] == 0 && blockIdx[end] == length(x0)
                        nnz = length(jac_g_colind) > 0 ? length(jac_g_row) : -1
                        
                        @assert length(vgroup_types) == length(vgroup_ind)
                        
                        new(Cint(length(lb_var)), Cint(length(lb_con)), Cint(nnz), [Cint(x) for x in blockIdx], 
                            vblocks, condenser, vgroup_types, vgroup_ind .|> x -> x .|> Cint,
                            [Cdouble(x) for x in lb_var], [Cdouble(x) for x in ub_var], [Cdouble(x) for x in lb_con], [Cdouble(x) for x in ub_con], Cdouble(lb_obj), Cdouble(ub_obj),
                            f, g, grad_f, jac_g, last_hessBlock, hess,
                            # constrVioReducer,
                            jac_g_nz, [Cint(x) for x in jac_g_row], [Cint(x) for x in jac_g_colind], 
                            x0, lambda0
                            )
    end
end

function make_sparse!(B_prob::Problem, nnz::Integer, jac_nz::Function, jac_row::Vector{T}, jac_col::Vector{T}) where T <: Integer
    B_prob.jac_g_nz = jac_nz
    B_prob.jac_g_row = [Cint(x) for x in jac_row]
    B_prob.jac_g_colind = [Cint(x) for x in jac_col]
    B_prob.nnz = Cint(nnz)
end
