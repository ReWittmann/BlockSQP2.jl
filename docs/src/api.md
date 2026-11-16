# API reference

The short names `Problem`, `Options`, `Stats`, `Solver`, and `Optimizer` are
available with the `BlockSQP2.` qualifier. Exported `BlockSQP2Problem`,
`BlockSQP2Options`, `BlockSQP2Stats`, `BlockSQP2Solver`, and `BlockSQP2Optimizer`
are aliases for the same types.

## Solver and condensing

```@autodocs
Modules = [BlockSQP2]
```

## Result codes

```@docs
BlockSQP2.SQPresults
```

## NLP layouts

Use `using BlockSQP2.NLPlayouts` to access layout helpers, ComponentArrays, and the
block kind aliases `nlpBlock`, `nlpVariables`, `nlpConstraints`, `nlpHess`,
`nlpMatching`, `nlpMatchings`, `nlpMultipleShootingDF`, `nlpMSfree`, and
`nlpMSdependent`. The unprefixed abstract types below describe those same kinds.

```@autodocs
Modules = [BlockSQP2.NLPlayouts]
```

```@index
```
