-- SPDX-License-Identifier: MIT
-- Copyright (c) 2026-present K. S. Ernest (iFire) Lee

import Init.Data.FloatArray
import Init.Data.Array

-- ============================================================================
-- 1. DOMAIN & DATA STRUCTURES
--
-- All spatial coordinates are in MICROMETRES (Int).
-- All velocities are in μm/tick (Int, abs, non-negative after FFI conversion).
-- All accelerations are in μm/tick² (Int, abs, non-negative after FFI conversion).
-- tickHz is a parse-time parameter; the AST never sees it after parseLeaf.
-- All costs are in μm² (Int).
--
-- Converting at the FFI boundary (parseLeaf) keeps the internal algorithm
-- entirely in ℤ, making every comparison decidable and every property
-- provable without Float lemmas or sorry.
-- ============================================================================

-- #snippet domain-aliases
/-- Identifier for an equivalence class in the spatial E-graph. -/
abbrev EClassId := Nat
/-- Identifier for a partition node in the spatial E-graph. -/
abbrev ENodeId  := Nat
-- #end domain-aliases

-- Axis-aligned bounding box; all coordinates in micrometres.
-- #snippet BoundingBox
/-- Axis-aligned bounding box with all coordinates in integer micrometres. -/
structure BoundingBox where
  minX : Int
  maxX : Int
  minY : Int
  maxY : Int
  minZ : Int
  maxZ : Int
  deriving Inhabited, Repr, DecidableEq
-- #end BoundingBox

/-- Integer-valued 3D vector (micrometres). Used by Plane and by segment/convex queries. -/
structure Vec3 where
  x : Int
  y : Int
  z : Int
  deriving Inhabited, Repr, DecidableEq

/-- Oriented half-space `{p : normal · p + d ≥ 0}`. The *kept* side is where
    `normal · p + d ≥ 0`; a point with `normal · p + d < 0` is rejected.
    Matches DynamicBVH's `Plane::is_point_over` convention but inverted
    so `is_point_over p = ¬ keeps p`. -/
structure Plane where
  normal : Vec3
  d      : Int
  deriving Inhabited, Repr, DecidableEq

-- Physics leaf: one entity, all quantities in μm / μm/tick / μm/tick².
-- timeOffset: ticks from "now" when this ghost's window starts (0 = current).
-- duration:   ticks this ghost covers (0 = use the BVH's global δstar).
-- Ghost leaves created by recursiveGhostSplit set both fields; parseLeaf leaves both 0.
-- #snippet LeafData
/-- Physics leaf representing one entity: bounds in μm, velocity in μm/tick, acceleration in μm/tick².
    Ghost leaves set timeOffset and duration; regular parseLeaf leaves leave both at 0. -/
structure LeafData where
  entities     : Nat
  bounds       : BoundingBox
  velocity     : Array Int   -- [|Vx|, |Vy|, |Vz|]  in μm/tick  (abs; = |v_μm_per_s| / tickHz)
  acceleration : Array Int   -- [|Ax|, |Ay|, |Az|]  in μm/tick² (abs; = ⌈½|a_μm_per_s2| / tickHz²⌉)
  timeOffset   : Nat := 0    -- ticks from "now" when this ghost window starts
  duration     : Nat := 0    -- ticks this ghost covers (0 = global δstar)
  deriving Inhabited, Repr
-- #end LeafData

-- #snippet PartitionNode
-- 3D partition node vocabulary.
-- AV1 SPACE-FILLING INVARIANT:
-- Every emitted split node stores `parent : BoundingBox` — the union of entity-tight
-- child bounds at the time the node was created.  Child geometric bounds are derived
-- as exact midpoint halves/quarters/eighths of `parent`.  The stored `parent` is the
-- proof witness used by the space-filling theorems in Section 8d: for any axis-aligned
-- midpoint split, the child regions tile `parent` exactly (manifold + space-filling).
-- Entity-tight bounds remain on leaf EClasses for SAH quality; `parent` enables proofs.
/-- 3D BVH partition node vocabulary covering 2-, 3-, 4-, 8-way splits and leaf.
    Every split node stores a `parent` BoundingBox (the octree cell) as the space-filling proof witness. -/
inductive PartitionNode where
  | none_split (data : LeafData)
  | horz       (parent : BoundingBox) (top bot : EClassId)      -- Y-axis 2-way
  | vert       (parent : BoundingBox) (left right : EClassId)   -- X-axis 2-way
  | depth      (parent : BoundingBox) (front back : EClassId)   -- Z-axis 2-way
  | horz_a     (parent : BoundingBox) (tl tr bot : EClassId)    -- XY T-shape: top split X, bottom full
  | vert_b     (parent : BoundingBox) (left tr br : EClassId)   -- XY T-shape: left full, right split Y
  | xz_a       (parent : BoundingBox) (fl fr back : EClassId)   -- XZ T-shape: front split X, back full
  | horz_4     (parent : BoundingBox) (s1 s2 s3 s4 : EClassId)  -- 4 Y-strips
  | vert_4     (parent : BoundingBox) (s1 s2 s3 s4 : EClassId)  -- 4 X-strips
  | depth_4    (parent : BoundingBox) (s1 s2 s3 s4 : EClassId)  -- 4 Z-strips
  | oct        (parent : BoundingBox) (s1 s2 s3 s4 s5 s6 s7 s8 : EClassId)
  deriving Inhabited, Repr
-- #end PartitionNode

-- ============================================================================
-- 2. E-GRAPH STORAGE
-- ============================================================================

-- #snippet EClass
/-- An equivalence class in the spatial E-graph: tracks the best (lowest-cost) partition node,
    entity-tight bounds for SAH, and Morton-code range for cell reconstruction. -/
structure EClass where
  id        : EClassId
  nodes     : Array ENodeId
  minCost   : Int              -- SAH cost in μm²
  bestNode  : Option ENodeId
  bounds    : BoundingBox      -- entity-tight bounds for SAH union computation
  firstCode : Nat              -- Morton code of leftmost (first) entity in this class
  lastCode  : Nat              -- Morton code of rightmost (last) entity in this class
  deriving Inhabited
-- #end EClass

-- #snippet SpatialEGraph
/-- The full spatial E-graph: flat arrays of nodes and classes, the scene AABB,
    the optimal prediction window δ*, and the root class after saturation. -/
structure SpatialEGraph where
  nodes        : Array PartitionNode
  classes      : Array EClass
  rootId       : Option EClassId    -- set by applyRewrites after saturation
  scene        : BoundingBox        -- full scene AABB; used for Morton cell reconstruction
  optimalDelta : Nat                -- auto-computed optimal prediction window (ticks); 1 = rebuild every tick
  deriving Inhabited
-- #end SpatialEGraph

-- ============================================================================
-- 2b. BOUNDING BOX UTILITIES
-- ============================================================================

-- #snippet surfaceArea
/-- Surface area 2(wh + hd + wd) of a bounding box in μm². Used as the SAH traversal cost proxy. -/
def surfaceArea (b : BoundingBox) : Int :=
  let w := b.maxX - b.minX
  let h := b.maxY - b.minY
  let d := b.maxZ - b.minZ
  2 * (w * h + h * d + w * d)
-- #end surfaceArea

/-- Tight axis-aligned union of two bounding boxes; used for SAH parent computation. -/
def unionBounds (a b : BoundingBox) : BoundingBox :=
  { minX := min a.minX b.minX, maxX := max a.maxX b.maxX,
    minY := min a.minY b.minY, maxY := max a.maxY b.maxY,
    minZ := min a.minZ b.minZ, maxZ := max a.maxZ b.maxZ }

/-- BVH traversal coefficient (dimensionless): 1 means traversal of a parent costs 1 μm² leaf-equivalent. -/
def bvhTraversalCost : Int := 1

/-- The single source of truth for simulation tick rate.
    20 Hz = current default. 10 Hz = hard floor (VRChat IK sync rate; also the
    ≤100 ms Long-Latency-Reflex bound).
    Values below 10 Hz break the mocap-freshness guarantee — do not lower.
    Change this ONE value to retarget the entire pipeline. -/
def simTickHz : Nat := 20

/-- Interest radius: 5 m = 5,000,000 μm. Single source of truth for ghost snap δ
    computation (CodeGen) and adversarial sim (Sim.lean). Internal unit: μm. -/
def interestRadius : Nat := 5000000

/-- Physical velocity ceiling: human sprint ≈ 10 m/s.
    Per-tick value depends on simTickHz: 10 m/s × 1,000,000 μm / tickHz.
    Server-enforced maximum for ghost bound clamping (C1/G13 fix). -/
def vMaxPhysicalAt (tickHz : Nat) : Nat := 10 * 1000000 / (max tickHz 1)

/-- vMaxPhysical at the configured simTickHz. -/
def vMaxPhysical : Nat := vMaxPhysicalAt simTickHz

/-- Hysteresis threshold: simulation steps an entity must remain in the new zone before STAGING begins.
    4 seconds expressed as step count (`simTickHz * 4`). Single source of truth —
    emitted as `pbvh_hysteresis_threshold(hz)` in the generated C/Rust. -/
def hysteresisThreshold : Nat := simTickHz * 4

-- ============================================================================
-- SPACE-FILLING CURVE UTILITIES (for O(N+k) broadphase)
-- ============================================================================

/-- Count leading zeros for a 30-bit space-filling curve code (result ∈ [0, 30]).
    Curve-agnostic: the width is the code's, not the curve's. -/
def clz30 (x : Nat) : Nat :=
  if x == 0 then 30 else 29 - Nat.log2 x

-- ── 3D Morton curve (Z-order) ─────────────────────────────────────────
--
-- This was a Hilbert curve until 2026-08-12 and is now Morton. The reason is not that
-- Hilbert is worse -- measured on this workspace's own data it is better on the metric it
-- is chosen for, 568 disjoint query ranges against Morton's 868, with a worst-case run
-- extent bounded at 2*L^(1/3)-1 where Morton has no bound at all.
--
-- It was written wrong twice. Here, and again by hand in C in the Godot engine's
-- `core/math/predictive_bvh_adapter.h`, with the same two mistakes both times: Skilling's
-- inner loop missing its `i = 0` step, and the transpose emitting z,y,x where X[0]=x
-- belongs. Measured, both scored 3583 of 4095 consecutive codes NOT face-adjacent -- 87.5%,
-- which is worse locality than the Morton they were chosen over, at five times the cost.
--
-- Both were clean bijections. Both round-tripped. The C one carried a CRASH_COND round-trip
-- witness that fired on every call and never once failed. A paper was cited. None of that
-- separates a Hilbert curve from any other bijection, and the defining property -- that
-- consecutive codes are adjacent cells -- was asserted nowhere until it was too late.
--
-- Morton cannot fail that way, and that is the whole argument. It is a bit permutation, so
-- the octree-prefix property the zone partitioning actually depends on is its definition
-- rather than a consequence of getting a rotation state machine right. There is no state to
-- get backwards. For prefix-partitioned assignment the two curves measured EXACTLY equal on
-- seam cost at every power-of-two zone count, so the locality Hilbert genuinely wins is not
-- being spent where this code spends it.

/-- Spread the low 10 bits of `n` so that bit i lands at position 3i. -/
def part1by2 (n : Nat) : Nat :=
  let n := n &&& 0x3ff
  let n := (n ||| (n <<< 16)) &&& 0x030000FF
  let n := (n ||| (n <<< 8)) &&& 0x0300F00F
  let n := (n ||| (n <<< 4)) &&& 0x030C30C3
  (n ||| (n <<< 2)) &&& 0x09249249

/-- Inverse of `part1by2`: gather every third bit back down. -/
def compact1by2 (n : Nat) : Nat :=
  let n := n &&& 0x09249249
  let n := (n ||| (n >>> 2)) &&& 0x030C30C3
  let n := (n ||| (n >>> 4)) &&& 0x0300F00F
  let n := (n ||| (n >>> 8)) &&& 0x030000FF
  (n ||| (n >>> 16)) &&& 0x000003FF

/-- A 30-bit 3D Morton index from three 10-bit coordinates.

    `x` leads each triple. Two implementations disagreeing about that produce different
    curves while both round-tripping perfectly, which is exactly how the previous defect
    stayed hidden -- so it is stated here and gated below rather than left to convention. -/
def morton3D (x y z : Nat) : Nat :=
  (part1by2 x <<< 2) ||| (part1by2 y <<< 1) ||| part1by2 z

/-- 30-bit code back to three 10-bit coordinates. -/
def morton3DInverse (m : Nat) : Nat × Nat × Nat :=
  (compact1by2 (m >>> 2), compact1by2 (m >>> 1), compact1by2 m)

/-- Transitional alias. `hilbert3D` IS `morton3D` and has been since 2026-08-12.

    The name is kept only so the ~240 references across `lean-spatial-oracle`,
    `lean-rebac-core` and the Godot module do not all have to move in one commit. Every
    use site gets a deprecation warning, so the compiler carries the remaining rename
    rather than a TODO nobody reads.

    This is a lying name on purpose and for a bounded time, which is a different thing
    from what went wrong here before: that was a function claiming to be Hilbert and
    being a BROKEN Hilbert, with nothing saying so. Delete this alias once the callers
    have moved. -/
@[deprecated morton3D (since := "2026-08-12")]
def hilbert3D (x y z : Nat) : Nat := morton3D x y z

-- ── The gate ───────────────────────────────────────────────────────────
--
-- These replace the contiguity gate the Hilbert version carried. Morton would FAIL that one
-- by design -- half its consecutive codes jump -- so keeping it would have been a test that
-- passes for the wrong curve, which is the failure being fixed rather than repeated.
--
-- What is gated instead is what the partitioning actually relies on: the code is a bijection,
-- and the top 3d bits name the octree cell of side 2^(10-d).

/-- Every cell of an `n`-cube, as `(x, y, z)`. -/
def cubeCells (n : Nat) : List (Nat × Nat × Nat) :=
  (List.range n).flatMap fun x => (List.range n).flatMap fun y => (List.range n).map fun z => (x, y, z)

/-- Codes that do not agree with their own octree cell at depth `d`. This is THE property:
    a prefix of `3*d` bits must name the cell of side `2^(10-d)`, because that is what a
    zone span is. -/
def mortonPrefixDefects (n d : Nat) : Nat :=
  let shift := if d ≤ 10 then 10 - d else 0
  ((cubeCells n).filter (fun c =>
    let code := morton3D c.1 c.2.1 c.2.2
    let origin := morton3D ((c.1 >>> shift) <<< shift) ((c.2.1 >>> shift) <<< shift)
      ((c.2.2 >>> shift) <<< shift)
    ¬ (code >>> (3 * shift) = origin >>> (3 * shift)))).length

theorem morton_roundtrips :
    ((cubeCells 16).filter (fun c =>
      ¬ (morton3DInverse (morton3D c.1 c.2.1 c.2.2) = (c.1, c.2.1, c.2.2)))).length = 0 := by
  native_decide

theorem morton_is_a_bijection :
    ((cubeCells 16).map (fun c => morton3D c.1 c.2.1 c.2.2)).eraseDups.length = 4096 := by
  native_decide

/-- THE GATE. The prefix is the octree cell, at every depth. -/
theorem morton_prefix_is_the_octree_cell :
    mortonPrefixDefects 16 6 = 0 ∧ mortonPrefixDefects 16 8 = 0 ∧
      mortonPrefixDefects 16 10 = 0 := by native_decide

