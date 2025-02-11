function lower_to_full!(arr1::Array{Float64, 1}, arr2::Array{Float64, 1}, n::Int32)
    for i = 1:n
        for j = 0:(i-1)
            arr2[i + j*n] = arr1[i + j*n - Int64(j*(j+1)//2)]
        end
    end

    for i = 1:n
        for j = i:(n-1)
            arr2[i+j*n] = arr1[(j+1) + (i-1)*n - Int64(i*(i-1)//2)]
        end
    end
end

function full_to_lower!(arr1::Array{Float64, 1}, arr2::Array{Float64, 1}, n::Int32)
    for i = 1:n
        for j = 0:(i-1)
            arr2[i + j*n - Int64(j*(j+1)//2)] = arr1[i + j*n]
        end
    end
end

function make_sparse(B_prob::blockSQPProblem, nnz::Int32, jac_nz::Function, jac_row::Array{Int32, 1}, jac_col::Array{Int32, 1})
    B_prob.jac_g_nz = jac_nz
    B_prob.jac_g_row = jac_row
    B_prob.jac_g_col = jac_col
    B_prob.nnz = nnz
end


function reduceConstrVio(Prob::Ptr{Nothing}, xi::CxxPtr{Float64}, info::CxxPtr{Int32})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    if Jprob.continuity_restoration == fnothing
        info[] = Int32(1)
    else
        xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
        xi_arr[:] = Jprob.continuity_restoration(xi_arr)
    end
    return
end