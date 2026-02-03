"""
    System of NLP variables resulting from a multiple shooting,
    with the variables being ordered as DF [Dependent,Free] in each
    shooting interval. All direct children must be of blocktype Hess,
    and have subblocks indicating dependent and free sections.
    
    MultipleShootingDF
    |...Hess
        |...MSfree
    |...Hess
        |...MSdependent
        |...MSfree
    |...Hess
     .
     .
     .
    |...Hess
        |...MSdependent
        |...MSfree
    
    The dependent a free sections may have further sublayouts, marking
    e.g. controls and parameters.
    In addition, it must have an attribute :matchings referring to a 
    MultipleShootingMatchings constraint block indicating the state matching
    constraints of the system.
"""
abstract type MultipleShootingDF <: Variables end

abstract type MSfree <: Variables end
abstract type MSdependent <: Variables end



# abstract type ParameterMatching <: Matching end
# abstract type ParameterMatchings <: Matchings end
# abstract type ControlMatchings <: Matchings end


function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: MultipleShootingDF
    @assert :matchings in keys(kwargs) && blocktypeof(kwargs[:matchings]) <: Matchings
    return BlockDescriptor{Btype{MultipleShootingDF}}(args...; kwargs...)
end


function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: MultipleShootingDF
    MP = tagmap(struc)
    axsubBD(AX, ind) = let __MP = MP
        axsubkeys(AX, ind) .|> x->__MP[x]
    end
    
    #a) Direct subblocks must be Hessian blocks
    Hchildren = axsubBD(struc.vLayout, BD)
    @assert all(isa.(blocktypeof.(Hchildren), Type{Hess})) "All direct sub-BlockDescriptors of a MultipleShootingSystemSC must have BlockType \"Hess\""
    
    #b) One matching block for each shooting stage
    MTC = BD.attr[:matchings]
    Mchildren = axsubBD(struc.cLayout, MTC)
    @assert (length(Mchildren) == length(Hchildren) - 1) "Number matching conditions does not match Hessian block structure"
    
    #c) First Hessian block only has one subblock of type MSfree
    DFchildren = axsubBD(struc.vLayout, first(Hchildren))
    # @assert (length(axsubkeys(struc.vLayout, first(Hchildren))) == 1) "First Hessian block of a MultipleShootingSystemSC must have exactly one subblock"
    @assert (length(DFchildren) == 1 && blocktypeof(first(DFchildren)) == MSfree) "First Hessian block of a MultipleShootingSystemSC must have exactly one subblock"

    
    #d) Hessian blocks except the first and the last must have one control/parameter C subblock and one state S subblock
    for i in eachindex(Hchildren)[2:end-1]
        # subtags = axsubkeys(struc.vLayout, Hchildren[i])
        DFchildren = axsubBD(struc.vLayout, Hchildren[i])
        # @assert length(subtags) == 2 "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (free section and dependent section)"
        @assert (length(DFchildren) == 2 && blocktypeof(first(DFchildren)) == MSdependent && blocktypeof(last(DFchildren)) == MSfree) "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (MSfree and MSdependent)"

        # @assert (length(axsubrange(struc.vLayout, MP[first(subtags)])) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
        @assert (length(axsubrange(struc.vLayout, first(DFchildren))) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
    end
    
    #e) Last Hessian block may have one dependent or one dependent and one free subblock.
    i = lastindex(Hchildren)
    DFchildren = axsubBD(struc.vLayout, Hchildren[i])
    @assert (1 <= length(DFchildren) <= 2 && blocktypeof(first(DFchildren)) == MSdependent && (length(DFchildren) == 1 || blocktypeof(last(DFchildren)) == MSfree)) "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (MSfree and MSdependent)"
    @assert (length(axsubrange(struc.vLayout, first(DFchildren))) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
    
    return nothing
end

# function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: Matchings
#     MP = tagmap(struc)
#     axsubBD(AX, ind) = let __MP = MP
#         axsubkeys(AX, ind) .|> x->__MP[x]
#     end
    
#     for subBD in axsubBD(struc.cLayout, BD)
#         @assert (blocktypeof(subBD) == Matching) "Direct subblocks of a Matchings must be of blocktype Matching"
#     end
# end
