# Guide: Interpreting Assembly Solver Inconsistent Constraint Diagnostics

## Overview

When the assembly solver cannot find a valid configuration, it produces an `INCONSISTENT_CONSTRAINTS` diagnostic. This means the constraints are **geometrically impossible** to satisfy simultaneously—not merely redundant or under-specified.

This guide is designed for an AI agent that will analyze the diagnostic, interpret joint visualizations, and make iterative corrections to the part and assembly JSON files to resolve the issues.

**Important**: Ground joints (Locked type connecting a part to the assembly origin) are automatically filtered from the diagnostic. They define the reference frame and cannot be "inconsistent" by themselves. If a ground joint would appear violated, it means other joints are pulling the part away from where it was grounded—the conflict is with those other joints.

---

## What You Receive

When analyzing an inconsistent assembly, you receive:

1. **Diagnostic JSON**: Parsed from the solver's YAML output, containing:
   - `affected_part`: The part appearing in most violations
   - `total_violation`: Sum of all constraint violations (lower = better)
   - `joints`: Array of joints with violations and LCS positions
   - `part_constraints`: Any part-level constraint violations

2. **Joint Visualization Images**: For each joint in the diagnostic:
   - Shows both connected parts and their LCS orientations
   - LCS axes displayed as colored lines (Red=X, Green=Y, Blue=Z)
   - Helps visualize misalignments that correspond to violation values

3. **Part JSON Files**: Part definitions you can modify
4. **Assembly JSON File**: Joint specifications you can modify

---

## What You Can Modify

### Part JSON: `dimensions`

```json
{
  "part_name": "PostFrontLeft",
  "dimensions": {
    "length": 80,   // EDITABLE
    "width": 80,    // EDITABLE
    "height": 600   // EDITABLE
  },
  "material": "oak",
  "theta_resolution": 90.0
}
```

**When to modify**: Part too short/long to span required distance.

### Assembly JSON: `refs` (Keypoint Selection)

```json
{
  "refs": ["CenterOfGravity_X40_Y40_Z0", "CenterOfGravity_X0_Y0_Z300"]
}
```

Keypoint IDs follow the pattern: `<Type>_X<x>_Y<y>_Z<z>`
- Examples: `CenterOfGravity_X50_Y25_Z100`, `Vertex_X0_Y0_Z0`, `ArcCenter_X100_Y50_Z50`

**When to modify**: Wrong attachment point selected, or need different surface normal direction.

### Assembly JSON: `theta` (Rotation Angle)

```json
{
  "theta": 0  // EDITABLE - rotation in degrees
}
```

- First part uses `+` LCS (Z pointing outward from surface)
- Second part uses `-_theta` LCS (anti-parallel to first, rotated by theta)
- `theta` must be a multiple of the part's `theta_resolution` (default: 90)
- Common values: 0, 90, 180, 270

**When to modify**: Orientation mismatch between parts.

### Assembly JSON: Add/Remove Joints

You can add new joints or remove existing ones from the `joints` array.

**When to modify**: Over-constrained (too many joints) or under-constrained assembly.

---

## Diagnostic JSON Structure

```json
{
  "affected_part": "PartName",
  "total_violation": 0.123456,
  "inconsistent_equation_count": 3,
  "nqsu": 14,
  "joints": [
    {
      "type": "FixedJoint",
      "part_i": "Part1",
      "part_j": "Part2",
      "lcs_i": {
        "name": "LCS_Origin",
        "position_on_part": [0.0, 0.0, 0.0],
        "world_position": [1.0, 2.0, 3.0],
        "world_quaternion": [0.0, 0.0, 0.0, 1.0]
      },
      "lcs_j": {
        "name": "LCS_Origin",
        "position_on_part": [0.0, 0.0, 0.0],
        "world_position": [1.1, 2.0, 3.0],
        "world_quaternion": [0.0, 0.0, 0.0, 1.0]
      },
      "relative_angle_degrees": 0.5,
      "relative_quaternion": [0.0, 0.0, 0.004, 1.0],
      "inconsistent_constraints": [
        {"type": "TranslationConstraintIJ directioni", "violation": 0.1}
      ]
    }
  ],
  "part_constraints": [
    {"part": "Part1", "type": "EulerParameterConstraint", "violation": 0.0001}
  ]
}
```

---

## Field Descriptions

### Top-Level Fields

| Field | Description |
|-------|-------------|
| `affected_part` | Part appearing most frequently in violated constraints. Often (but not always) the misconfigured part. |
| `total_violation` | Sum of absolute violations. **Use this to measure improvement between iterations.** |
| `inconsistent_equation_count` | Number of constraint equations that could not be satisfied. |
| `nqsu` | Number of degrees of freedom variables. |

### Joint Fields

| Field | Description |
|-------|-------------|
| `type` | Joint class (FixedJoint, RevoluteJoint, etc.) |
| `part_i`, `part_j` | The two parts connected by this joint |
| `relative_angle_degrees` | Angular misalignment between the two LCS frames (0° = perfectly aligned) |
| `relative_quaternion` | Relative rotation as quaternion `conj(qI) * qJ` |

### LCS Fields

| Field | Description |
|-------|-------------|
| `name` | LCS identifier |
| `position_on_part` | LCS origin relative to part's origin (mm) |
| `world_position` | LCS origin in global coordinates (mm) |
| `world_quaternion` | LCS orientation as quaternion [x, y, z, w] |

---

## Constraint Types and Their Meanings

### Position Constraints (violations in mm)

| Type | Meaning | Violation Interpretation |
|------|---------|-------------------------|
| `AtPointConstraintIJx/y/z` | Points must coincide on X/Y/Z axis | Distance apart on that axis |
| `TranslationConstraintIJ directioni/j/k` | Translation component constraint | Distance error on that direction |
| `DispCompIecJecO` | Displacement component constraint | Distance error |

### Orientation Constraints (violations are direction cosines, range -1 to +1)

| Type | Meaning | Violation Interpretation |
|------|---------|-------------------------|
| `DirectionCosineConstraintIxJx` | X-axes must align | 0 = perpendicular, ±1 = aligned |
| `DirectionCosineConstraintIxJy` | X of I perpendicular to Y of J | Should be 0; non-zero = not perpendicular |

**Converting direction cosine violations to angles:**
- `violation = 0` → axes are correctly perpendicular (90°)
- `violation = ±0.45` → axes are at ~63° instead of 90° (27° error)
- `violation = ±1.0` → axes are parallel instead of perpendicular (90° error)

---

## Using `relative_angle_degrees` for Quick Diagnosis

| Value | Interpretation |
|-------|----------------|
| 0° - 1° | Frames nearly aligned; issue is likely positional |
| 1° - 10° | Small angular misalignment; try adjusting `theta` |
| 10° - 45° | Significant misalignment; wrong `refs` or `theta` |
| 45° - 90° | Major misalignment; part may be rotated 90° wrong |
| > 90° | Severe misalignment; try `theta` ± 180° |

---

## Interpreting Joint Visualizations

The visualization images show each joint with both connected parts and their LCS orientations.

### Reading LCS Axes

- **Red axis**: X direction
- **Green axis**: Y direction
- **Blue axis**: Z direction (normal to surface)

### What to Look For

**Correct alignment (Fixed joint)**:
- Blue (Z) axes of both LCS should point toward each other (anti-parallel)
- The LCS origins should be at the same world position
- X and Y axes should align based on the `theta` value

**Position issues**:
- Visible gap between LCS origins indicates position violation
- Large gaps suggest wrong `refs` selection or incorrect part dimensions

**Orientation issues**:
- Blue axes not anti-parallel indicates orientation mismatch
- Try changing `theta` by 90° or 180°
- If axes point in completely unexpected directions, wrong `refs` keypoint

### Correlating Visuals with Violations

| Visual Observation | Likely Cause | Diagnostic Indicator |
|-------------------|--------------|---------------------|
| Gap between LCS origins | Wrong keypoint or part dimension | High `TranslationConstraint` violations |
| Z-axes not anti-parallel | Wrong theta value | High `DirectionCosine` violations |
| Axes at ~90° to expected | theta off by 90° | `relative_angle_degrees` ≈ 90 |
| Parts completely misoriented | Wrong refs keypoint | Large `relative_angle_degrees` |

---

## Diagnostic Reasoning Process

### Step 1: Analyze Joint Visualizations

Examine the visualization images for each joint in the diagnostic:
- Are the LCS origins at the same location?
- Are the Z-axes pointing toward each other?
- Does the alignment match what you'd expect for the intended assembly?

### Step 2: Identify the Conflict Pattern

**Single joint with violations:**
- The joint's two LCS frames have incompatible positions or orientations
- Check if one part needs different `refs` or `theta`

**Multiple joints with violations involving a common part:**
- That part is being "pulled" in incompatible directions
- The part cannot satisfy all joints simultaneously
- Consider: Is this part's geometry wrong? Are the refs correct?

**Chain/loop of joints:**
- Closed kinematic loop with incompatible geometry
- Errors accumulate around the loop
- May need to adjust intermediate joints' theta values

### Step 3: Analyze World Positions

Compare `world_position` of `lcs_i` and `lcs_j` within each joint:

- **Positions match**: Issue is purely rotational → adjust `theta`
- **Positions differ**: Position constraint violated → change `refs` or `dimensions`

**Large position differences** (hundreds of mm) suggest wrong keypoint selection or grossly incorrect part dimensions.

### Step 4: Analyze Violation Magnitudes

| Violation Size | Interpretation |
|---------------|----------------|
| < 0.01 | Near-tolerance; may resolve with small theta adjustment |
| 0.01 - 0.5 | Significant conflict; likely wrong refs or theta |
| > 0.5 | Major conflict; possibly wrong keypoint face or part dimension |

### Step 5: Check `relative_angle_degrees`

Quick summary of angular misalignment per joint:
- High values → orientation problem → change `theta`
- Low values with position violations → translation problem → change `refs` or `dimensions`

---

## Common Failure Patterns

### Pattern 1: Orientation Mismatch in Chain
```
Joint A: Part1 → Part2, relative_angle_degrees = 26.5
Joint B: Part2 → Part3, relative_angle_degrees = 26.5
```
**Diagnosis**: Part2 (the common part) has an orientation that's ~26° off. It can't align with both Part1 and Part3.
**Fix**: Adjust theta on one of the joints by the error amount (if multiple of theta_resolution) or change refs.

### Pattern 2: Closed Loop Won't Close
```
Joint A: Part1 → Part2, relative_angle_degrees = 0.1
Joint B: Part2 → Part3, relative_angle_degrees = 0.1
Joint C: Part3 → Part1, relative_angle_degrees = 17.2
```
**Diagnosis**: Angular errors accumulate around the loop. When the loop tries to close, there's a ~17° mismatch.
**Fix**: Distribute theta adjustments across multiple joints, or fix the root cause (often wrong refs).

### Pattern 3: Positional Impossibility
```
Joint A: Part1 → Part2, TranslationConstraint violation = 50.0
Joint B: Part1 → Part3, TranslationConstraint violation = -50.0
```
**Diagnosis**: Part1 is being pulled 50mm in opposite directions.
**Fix**: Check if Part1's dimensions are correct, or if the refs keypoints are at wrong positions.

### Pattern 4: Part Too Short/Long
```
Joint A: Part1 → Part2, world_positions match
Joint B: Part1 → Part3, world_positions differ by [100, 0, 0]
```
**Diagnosis**: Part1's attachment points are 100mm apart in the model, but the assembly requires a different span.
**Fix**: Adjust Part1's `dimensions.length` (or width/height depending on orientation).

### Pattern 5: Anti-Parallel Configuration (Handled Automatically)
```
MbD: Protecting 2 essential constraint(s) from removal: 15 16
```
**Note**: The solver automatically handles anti-parallel (180°) configurations through iterative constraint protection. If you see these messages followed by successful convergence, the solver is working correctly.

---

## Making Corrective Changes

### For Position Violations

**Option 1: Change `refs` to different keypoint**
```json
// Before
"refs": ["CenterOfGravity_X40_Y40_Z0", "CenterOfGravity_X0_Y0_Z300"]

// After - select keypoint at different position
"refs": ["CenterOfGravity_X40_Y40_Z600", "CenterOfGravity_X0_Y0_Z300"]
```

**Option 2: Adjust part dimensions**
```json
// Before
"dimensions": {"length": 80, "width": 80, "height": 600}

// After - increase height to span required distance
"dimensions": {"length": 80, "width": 80, "height": 700}
```

### For Orientation Violations

**Primary fix: Adjust `theta`**
```json
// Before
"theta": 0

// After - try 90° increments
"theta": 90   // or 180, 270
```

**Alternative: Change refs to keypoint on different face**

If theta adjustments don't work, the keypoint may be on the wrong face. Select a keypoint with a different surface normal direction.

### For Common Part in Multiple Violations

The `affected_part` is often the culprit. Check:
1. Are all its refs keypoints on correct faces?
2. Are its dimensions correct?
3. Is it over-constrained (connected by too many joints)?

### For Over-Constrained Assembly

Remove redundant joints:
```json
// Before
"joints": [joint1, joint2, joint3, redundant_joint4]

// After
"joints": [joint1, joint2, joint3]
```

---

## Iterative Correction Workflow

Follow this process to systematically fix assembly issues:

```
1. ANALYZE
   - Review diagnostic JSON
   - Examine joint visualization images
   - Note total_violation value (baseline)

2. IDENTIFY
   - Find joint with highest violation
   - Determine if issue is positional or orientational
   - Check if affected_part appears in multiple violations

3. HYPOTHESIZE
   - Wrong refs? (position issue, wrong face)
   - Wrong theta? (orientation issue)
   - Wrong dimensions? (part too short/long)

4. CHANGE
   - Make SINGLE smallest change
   - Document what you changed and why

5. RE-RUN
   - Run solver again
   - Get new diagnostic

6. EVALUATE
   - Compare new total_violation to previous
   - If improved: continue with next issue
   - If worse: REVERT and try alternative hypothesis

7. REPEAT
   - Continue until total_violation = 0
   - Or until all constraints satisfied
```

**Critical**: Make only ONE change at a time. This allows you to identify which changes help and which hurt.

---

## Resolution Strategies

| Problem | Diagnostic Indicator | JSON Change |
|---------|---------------------|-------------|
| Wrong attachment point | High position violations | Change `refs` to keypoint at correct location |
| Parts at wrong angle | `relative_angle_degrees` ≈ 90 | Change `theta` by 90° |
| Parts facing wrong way | `relative_angle_degrees` ≈ 180 | Change `theta` by 180° |
| Part too short | Position violations, parts don't reach | Increase `dimensions.length/width/height` |
| Part too long | Position violations, parts overlap | Decrease `dimensions.length/width/height` |
| Over-constrained | Multiple redundant violations | Remove joint from `joints` array |
| Wrong face selected | Z-axes not anti-parallel | Change `refs` to keypoint on different face |

---

## Key Principles

1. **The solver reports symptoms, not root causes.** Multiple violations may stem from a single misconfigured part or joint.

2. **Start with the `affected_part`.** It appears in the most violations and is often the source of the problem.

3. **Use `relative_angle_degrees` for quick triage.** High values indicate orientation problems (fix with theta); low values suggest positional issues (fix with refs or dimensions).

4. **Make ONE change at a time.** This allows you to identify what helps and what hurts.

5. **Compare `total_violation` between iterations.** This is your progress metric. Lower = better.

6. **If a change makes things worse, revert immediately.** Try an alternative hypothesis.

7. **Prioritize highest-violation joints.** Fixing the biggest errors first often resolves smaller cascading issues.

8. **Look for common parts across multiple violations.** A part appearing in many violated joints is likely misconfigured.

9. **Orientation errors compound.** A 10° error at one joint becomes 20° by the time a chain reaches another joint.

10. **World positions reveal ground truth.** They show where parts actually ended up in the solver's best-effort state.

11. **Ground joints are filtered out.** If you see a common part in all violations, its ground joint is not the problem—the other joints are.

12. **Zero-violation constraints are filtered.** Only constraints that couldn't be satisfied appear in the diagnostic.
