# entities-lean-shared

Shared primitive types: the common vocabulary every `multiplayer-fabric` hexagon core builds on. Dependency-free (Lean core only).

> Split out of the [`lean-predictive-bvh`](https://github.com/v-sekai-multiplayer-fabric/lean-predictive-bvh) monorepo (now archived). Cross-cluster wiring is via Lake `require ... from git`.

## Dependencies

None — pure Lean (`Init`). Every other cluster requires this one.

## Build

```sh
lake build         # production gate: typecheck the Shared cluster
lake build Research  # research-tier (non-gating; may fail)
```

## Layout

- `Shared/` — the primitive types, in one Lean namespace

This repository has no `core/`, `ports/`, or `adapters/` directory. It holds a type vocabulary rather than a hexagon, so there is nothing on either side of a boundary to separate.
