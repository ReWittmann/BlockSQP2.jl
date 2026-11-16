# Building locally

The site uses [Documenter.jl](https://documenter.juliadocs.org/stable/man/guide/)
to render the manual and the package's Julia docstrings.

Use Julia 1.11 or later. Commands below run from the repository root.
This development checkout loads `BlockSQP2.jl/bin/libblockSQP2_jl` directly;
installing `blockSQP2_jll` alone does not replace that local library.
Build the native interface first using the dependencies and platform instructions
in the repository's main README:

```sh
mkdir -p BlockSQP2.jl/bin
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DJULIA_INTERFACE=ON
cmake --build build -j 4
```

Prepare the documentation environment and build HTML:

```sh
julia --project=BlockSQP2.jl/docs -e 'using Pkg; Pkg.develop(path="BlockSQP2.jl"); Pkg.instantiate()'
julia --project=BlockSQP2.jl/docs BlockSQP2.jl/docs/make.jl
```

Open `build/index.html` inside this directory. The site uses relative file URLs
and requires no server or deployment credentials. Documenter checks docstring
coverage and cross references, and executes the small getting started example.
The longer Lotka solve is a separate, headless check:

```sh
julia --project=BlockSQP2.jl/docs BlockSQP2.jl/docs/lotka.jl
```

`lotka.jl` adapts `../examples/example_lotka.jl`, omitting the Ipopt comparison,
variable-group scaling metadata, and interactive plotting. Keep its model and
layout in sync when editing the original example. Solver output uses a temporary
directory. Neither building nor running this script publishes anything.
