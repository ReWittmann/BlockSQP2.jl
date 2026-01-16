abstract type AbstractBlockDescriptor end

abstract type Block end

abstract type Variables <: Block end
abstract type Constraints <: Block end

abstract type Hess <: Variables end

abstract type Matching <: Constraints end

"""
Wrapper type for constructor selection
"""
abstract type Btype{T <: Block} <: Block end

mutable struct BlockDescriptor{T <: Block} <: AbstractBlockDescriptor
    tag::Symbol
    flags::Tuple
    attr::NamedTuple
    function BlockDescriptor{arg_T}(args...; tag = gensym(), kwargs...) where arg_T <: Block
        new{arg_T}(tag, (args...,), (; kwargs...))
    end
end

BlockDescriptor(args...; kwargs...) = BlockDescriptor{Block}(args...; kwargs...)

function BlockDescriptor{arg_BT}(args...; kwargs...) where arg_BT <: Btype{arg_T} where arg_T <: Block
    return BlockDescriptor{arg_T}(args...; kwargs...)
end

function blocktypeof(::BlockDescriptor{T}) where T
    return T
end


_tags(arg::AbstractVector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T} = arg .|> first .|> Base.Fix2(getfield, :tag)
