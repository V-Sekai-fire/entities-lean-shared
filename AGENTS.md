# AGENTS.md — lean-shared-core

Guidance for AI coding agents working in this repo.

## What this is

Shared primitive types — the common vocabulary every multiplayer-fabric hexagon core builds on. Dependency-free.

This is one hexagon cluster split out of the `lean-predictive-bvh` monorepo
(now archived). It follows the `core/ports/adapters` convention. Do not
implement algorithms elsewhere that contradict a proof here — if an
implementation differs from the Lean proof, trust the proof.

## Build

```sh
lake build           # production gate: `lake build Shared`
lake build Research  # research tier (non-gating; may fail against the pinned toolchain)
```

Lean toolchain version is pinned in `lean-toolchain`. Cross-cluster
dependencies are resolved via Lake `require ... from git` (see `lakefile.lean`);
run `lake update` to refresh `lake-manifest.json`.

## Conventions

- New algorithms land here first as Lean proofs, then get ported to Elixir/C++.
- Keep the production closure (the cluster aggregator) gating; park aspirational
  or toolchain-broken proofs under `Research.lean` (non-gating).
- Commit message style: sentence case, no `type(scope):` prefix.
  Example: `Prove O(1) refit bound for incremental BVH updates`
