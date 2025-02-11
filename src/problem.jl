mutable struct blockSQPProblem

    nVar::Int32
    nCon::Int32
    nnz::Int32
    blockIdx::Array{Int32, 1}

    lb_var::Array{Float64, 1}
    ub_var::Array{Float64, 1}
    lb_con::Array{Float64, 1}
    ub_con::Array{Float64, 1}
    objLo::Float64
    objUp::Float64

    f::Function
    g::Function
    grad_f::Function
    jac_g::Function
    jac_g_nz::Function
    continuity_restoration::Function
    last_hessBlock::Function
    hess::Function

    jac_g_row::Array{Int32, 1}
    jac_g_colind::Array{Int32, 1}

    x_start::Array{Float64, 1}
    lam_start::Array{Float64, 1}

    blockSQPProblem(nVar::Int32, nCon::Int32) = new(nVar, nCon, Int32(-1), Int32[0, nVar],
                    Float64[], Float64[], Float64[], Float64[], -Inf, Inf,
                    fnothing, fnothing, fnothing, fnothing, fnothing, fnothing, fnothing, fnothing,
                    Int32[], Int32[], Float64[], Float64[]
    )
end