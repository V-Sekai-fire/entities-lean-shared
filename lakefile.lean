-- SPDX-License-Identifier: MIT
-- Copyright (c) 2026-present K. S. Ernest (iFire) Lee

import Lake
open System Lake DSL

package «lean-shared-core» where

-- Shared primitive types: the common vocabulary every multiplayer-fabric
-- hexagon core builds on. Dependency-free (Lean core / Init only) so every
-- downstream cluster can require it without pulling Mathlib or any native FFI.
lean_lib Shared where
  roots := #[`Shared]
  globs := #[.one `Shared]
