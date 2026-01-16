module NLPstructures
using Reexport
@reexport using ComponentArrays
using DocStringExtensions

include("BlockDescriptor.jl")
export BlockDescriptor, blocktypeof

export Btype

nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching = (Btype{T} for T in (Block, Variables, Constraints, Hess, Matching))
export nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching

# _Block, _Variables, _Constraints, _Hess, _Matching = Block, Variables, Constraints, Hess, Matching
# export _Block, _Variables, _Constraints, _Hess, _Matching

#Note: The AbstractVector should be of type AbstractVector{TupleBD{S}}
TupleBD{S} = Tuple{BlockDescriptor, Union{AbstractVector, S}} where {S <: Integer}
export TupleBD


"""
$(TYPEDEF)
A generic type to store structure information of a nonlinear program,
e.g. Hessian block structure, dependency structure, etc.

# Fields
$(FIELDS)
"""
struct NLPstructure{VB, VL <: ComponentArrays.AbstractAxis, CB, CL <: ComponentArrays.AbstractAxis}
    "Collection of BlockDescriptor{Variables} describing blocks of variables"
    vBlocks::VB
    "A mapping of the BlockDescriptors to sections of the variables"
    vLayout::VL
    
    "Collection of all BlockDescriptor{Constraints} describing blocks of constraints"
    cBlocks::CB
    "See vLayout"
    cLayout::CL
end

tagmap(S::NLPstructure) = Dict((BD.tag => BD) for BD in union(S.vBlocks, S.cBlocks))

#Provide sanity checks for certain block types, e.g. MultipleShootingSystemSC
function is_valid(block::BlockDescriptor{B}, struc::NLPstructure) where B <: Block
    return true
end


export NLPstructure, tagmap


include("LayoutUtils.jl")
export axlength, axsubrange, to_NamedTuple, to_ComponentArray, to_UR, to_Axis, axsubindex, axsubkeys,
       simple_vBlocks, simple_cBlocks, hessBlocks, hessBlockSizes, hessBlockIndex, 
       has_parent, parent_of, has_parent_type


include("MultipleShootingSystemSC.jl")

nlpMultipleShootingSystemSC, nlpMultipleShootingMatchings = Btype{MultipleShootingSystemSC}, Btype{MultipleShootingMatchings}
export nlpMultipleShootingSystemSC, nlpMultipleShootingMatchings

# _MultipleShootingSystemSC, _MultipleShootingMatchings = MultipleShootingSystemSC, MultipleShootingMatchings
# export _MultipleShootingSystemSC, _MultipleShootingMatchings


end #module