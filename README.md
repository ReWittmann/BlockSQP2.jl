**BlockSQP2.jl** -- A Julia interface to blockSQP2, a nonlinear  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; programming solver based on [blockSQP](https://github.com/djanka2/blockSQP).  
Copyright (c) 2025 Reinhold Wittmann <reinhold.wittmann@ovgu.de> and Christoph Plate <christoph.plate@ovgu.de>

## Installation

This package requires the compiled binary of the C interface to [blockSQP2](https://github.com/ReWittmann/blockSQP2). One can either  
 - Fetch this package as a submodule of [blockSQP2](https://github.com/ReWittmann/blockSQP2), letting the build-process place the binary in BlockSQP2.jl/bin  
 - Depend on [blockSQP2_jll](https://github.com/JuliaRegistries/General/tree/master/jll/B/blockSQP2_jll)
  
The binary may link to, or include, compiled code subject to distinct licenses, see blockSQP2/README.md for a list of dependencies and their licenses.

Afterwards, the package can be managed through the Julia package manager as usual, see <https://docs.julialang.org/en/v1/stdlib/Pkg/>.

### Licensing
BlockSQP2.jl is published under the very permissive zlib license, see LICENSE.txt.
