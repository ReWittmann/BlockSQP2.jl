# module ComponentArraysExtension
using ComponentArrays
# using blockSQP
# @info "Loading ComponentArrays.jl extension for BlockSQP2..."

BlockSQP2.__lowerbounds(x::ComponentArray) = x[:]
BlockSQP2.__upperbounds(x::ComponentArray) = x[:]
BlockSQP2.__initial_values(x::ComponentArray) = x[:]

# end
