# entities-lean-shared

Shared Lean 4 primitive types that the multiplayer-fabric hexagon cores build on, with no dependency beyond the Lean core library.

## What it is for

It holds the geometric and spatial-partition types the cores have in common, with their helpers, in one namespace, so each core requires this library rather than defining its own. It is a vocabulary rather than a hexagon, so it has no core, ports or adapters.

## Build and run

```sh
lake build
```

## Licence

MIT; see `LICENSE`.
