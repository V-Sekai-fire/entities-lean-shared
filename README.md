# lean-shared-core

Shared primitive types: the common vocabulary every `multiplayer-fabric` hexagon core builds on. Dependency-free (Lean core only).

> Split out of the [`lean-predictive-bvh`](https://github.com/v-sekai-multiplayer-fabric/lean-predictive-bvh) monorepo (now archived). Each hexagon cluster is its own repo following the `core/ports/adapters` convention; cross-cluster wiring is via Lake `require ... from git`.

## Dependencies

None — pure Lean (`Init`). Every other cluster requires this one.

## Build

```sh
lake build         # production gate: typecheck the Shared cluster
lake build Research  # research-tier (non-gating; may fail)
```

## Hexagon layout

- `core/` — dependency-free domain logic + proofs
- `ports/` — narrow driving (source) / driven (sink) contracts
- `adapters/` — concrete I/O at the edges
