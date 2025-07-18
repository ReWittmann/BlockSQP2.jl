print(@__DIR__)

macro dlsym(lib, func)
    z = Ref{Ptr{Cvoid}}(C_NULL)
    quote
        let zlocal = $z[]
            if zlocal == C_NULL
                zlocal = Base.Libc.Libdl.dlsym($(esc(lib))::Ptr{Cvoid}, $(esc(func)))::Ptr{Cvoid}
                $z[] = zlocal
            end
            zlocal
        end
    end
end

libblockSQP = Base.Libc.Libdl.dlopen("../bin/libblockSQP_jl");

macro BSQP(func)
    fname = String(func)
    quote
        @dlsym(libblockSQP, $fname)
    end
end
