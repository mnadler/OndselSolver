# Guide: Interpreting Assembly Solver Inconsistent Constraint Diagnostics

## Overview

When the assembly solver cannot find a valid configuration, it logs an `INCONSISTENT_CONSTRAINTS` diagnostic to the console. This means the constraints are **geometrically impossible** to satisfy simultaneously—not merely redundant or under-specified.

**Important:** When constraints are inconsistent, the solver displays a **best-effort state** in FreeCAD's 3D view. This shows the closest configuration the solver could achieve before determining the constraints are impossible. This visual feedback helps identify which parts/joints are problematic.

## Message Structure

The diagnostic is wrapped with markers for programmatic parsing:

```yaml
---BEGIN:INCONSISTENT_CONSTRAINTS---
Constraints are geometrically inconsistent (no solution exists)
  affected_part: "<part_name>"           # Part appearing in most conflicting joints
  total_violation: <float>               # Sum of absolute violations
  inconsistent_equation_count: <int>     # Number of equations that couldn't be satisfied
  nqsu: <int>                            # Number of DOF variables (positions/orientations)
  joints:
    - name: "<joint_name>"
      type: "<joint_type>"               # e.g., "FixedJoint", "RevoluteJoint"
      part_i: "<part_name>"              # First part connected by this joint
      part_j: "<part_name>"              # Second part connected by this joint
      lcs_i:                             # Local Coordinate System on part_i
        name: "<lcs_name>"
        position_on_part: [x, y, z]      # LCS origin relative to part origin (mm)
        world_position: [x, y, z]        # LCS origin in global coordinates (mm)
        world_quaternion: [x, y, z, w]   # LCS orientation as quaternion [x, y, z, w]
      lcs_j:                             # Local Coordinate System on part_j
        name: "<lcs_name>"
        position_on_part: [x, y, z]
        world_position: [x, y, z]
        world_quaternion: [x, y, z, w]
      relative_angle_degrees: <float>    # Angle between lcs_i and lcs_j orientations
      relative_quaternion: [x, y, z, w]  # Relative rotation: conj(qI) * qJ
      inconsistent_constraints:
        - type: "<constraint_type>"
          violation: <float>             # How far from satisfied
  part_constraints:                      # Constraints on individual parts (e.g., Euler parameter)
    - part: "<part_name>"
      type: "<constraint_type>"
      equation_number: <int>
      violation: <float>
  unmatched_equations:                   # Equations not matched to known constraints
    - equation_number: <int>
      rhs_value: <float>
      hint: "<explanation>"
---END:INCONSISTENT_CONSTRAINTS---
```

## Field Descriptions

### Top-Level Fields

| Field | Description |
|-------|-------------|
| `affected_part` | The part that appears most frequently in violated constraints. Often (but not always) the misconfigured part. |
| `total_violation` | Sum of absolute values of all constraint violations. Larger = more severe conflict. |
| `inconsistent_equation_count` | Number of constraint equations that could not be satisfied. |
| `nqsu` | Number of degrees of freedom variables. Useful for understanding system size. |

### Joint Fields

| Field | Description |
|-------|-------------|
| `name` | Joint identifier from FreeCAD |
| `type` | Joint class (FixedJoint, RevoluteJoint, CylindricalJoint, etc.) |
| `part_i`, `part_j` | The two parts connected by this joint |
| `relative_angle_degrees` | Angular misalignment between the two LCS frames (0° = perfectly aligned) |
| `relative_quaternion` | Relative rotation as quaternion `conj(qI) * qJ` for precise orientation analysis |

### LCS (Local Coordinate System) Fields

| Field | Description |
|-------|-------------|
| `name` | LCS identifier (often "LCS_Origin") |
| `position_on_part` | LCS origin position relative to part's origin (mm) |
| `world_position` | LCS origin position in global/world coordinates (mm) |
| `world_quaternion` | LCS orientation in global coordinates as quaternion [x, y, z, w] |

### Part Constraints

These are constraints on individual parts, not joints:

| Field | Description |
|-------|-------------|
| `part` | Part name |
| `type` | Typically "EulerParameterConstraint" (ensures quaternion normalization) |
| `equation_number` | Internal equation index |
| `violation` | How far from satisfied |

### Unmatched Equations

Equations that couldn't be matched to known constraint types:

| Field | Description |
|-------|-------------|
| `equation_number` | Internal equation index |
| `rhs_value` | Right-hand side value at detection |
| `hint` | Explanation of what this might be (DOF equation vs unknown constraint) |

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
| `DirectionCosineConstraintIxJx` | X-axes must align | cosine of angle between axes (0 = perpendicular, ±1 = aligned) |
| `DirectionCosineConstraintIxJy` | X of I perpendicular to Y of J | Should be 0; non-zero = not perpendicular |
| `DirectionCosineConstraintIyJx` | Y of I perpendicular to X of J | Should be 0; non-zero = not perpendicular |
| `DirectionCosineConstraintIzJy` | Z of I perpendicular to Y of J | Should be 0; non-zero = not perpendicular |
| `DirectionCosineConstraintIzJx` | Z of I perpendicular to X of J | Should be 0; non-zero = not perpendicular |

**Converting direction cosine violations to angles:**
- `violation = 0` → axes are correctly perpendicular (90°)
- `violation = ±0.45` → axes are at arccos(0.45) ≈ 63° instead of 90° (27° error)
- `violation = ±1.0` → axes are parallel instead of perpendicular (90° error)

Formula: `angle_error = 90° - arccos(|violation|)`

### Part Constraints

| Type | Meaning | Violation Interpretation |
|------|---------|-------------------------|
| `EulerParameterConstraint` | Quaternion must be normalized (sum of squares = 1) | Deviation from unit length |

## Using `relative_angle_degrees` for Quick Diagnosis

The `relative_angle_degrees` field provides an immediate measure of how misaligned two LCS frames are:

| Value | Interpretation |
|-------|----------------|
| 0° - 1° | Frames nearly aligned; issue is likely positional |
| 1° - 10° | Small angular misalignment; minor adjustment needed |
| 10° - 45° | Significant misalignment; part orientation is wrong |
| 45° - 90° | Major misalignment; part may be rotated 90° wrong |
| > 90° | Severe misalignment; part may be flipped or facing wrong direction |

## Diagnostic Reasoning Process

### Step 1: Check the 3D View

The solver displays a **best-effort state** showing where parts ended up. Look for:
- Parts that are visibly misaligned
- Gaps between parts that should be connected
- Parts rotated to unexpected orientations

### Step 2: Identify the Conflict Pattern

**Single joint with violations:**
- The joint's two LCS frames have incompatible positions or orientations
- Check if one part is incorrectly positioned/rotated

**Multiple joints with violations involving a common part:**
- That part is being "pulled" in incompatible directions
- The part cannot satisfy all joints simultaneously

**Chain/loop of joints:**
- Closed kinematic loop with incompatible geometry
- The chain of parts doesn't close properly

### Step 3: Analyze World Positions

Compare `world_position` of `lcs_i` and `lcs_j` within each joint:

- **Positions match:** The attachment points are coincident; the issue is purely rotational
- **Positions differ:** Either:
  - A positional constraint is also violated, OR
  - The solver stopped before achieving positional convergence due to orientation impossibility

**Large position differences** (hundreds of mm) suggest a part is grossly misplaced or the wrong LCS is being used.

### Step 4: Analyze Violation Magnitudes

**Small violations (< 0.01):**
- Near-tolerance issues
- May indicate numerical precision problems or slightly over-constrained system

**Medium violations (0.01 - 0.5):**
- Significant geometric conflict
- Parts are oriented incorrectly by 1-30°

**Large violations (> 0.5):**
- Major geometric conflict
- Parts may be rotated ~45-90° from correct orientation
- Possibly wrong LCS selected or part inserted with wrong orientation

### Step 5: Check `relative_angle_degrees`

This gives a quick summary of angular misalignment per joint:
- High values indicate orientation problems
- Low values with position violations indicate translation problems

### Step 6: Look for Symmetry in Violations

**Violations with opposite signs (e.g., +0.45 and -0.45):**
- Suggests a part is rotated in one direction relative to where it should be
- The same angular error manifests as positive on one side, negative on the other

**All violations same sign:**
- Systematic error in one direction
- Likely a single part or LCS is incorrectly defined

## Common Failure Patterns

### Pattern 1: Orientation Mismatch in Chain
```
Joint A: Part1 → Part2, relative_angle_degrees = 26.5
Joint B: Part2 → Part3, relative_angle_degrees = 26.5
```
**Diagnosis:** Part2 (the common part) has an orientation that's ~26° off. It can't align with both Part1 and Part3.

### Pattern 2: Closed Loop Won't Close
```
Joint A: Part1 → Part2, relative_angle_degrees = 0.1
Joint B: Part2 → Part3, relative_angle_degrees = 0.1
Joint C: Part3 → Part1, relative_angle_degrees = 17.2
```
**Diagnosis:** The chain Part1→Part2→Part3→Part1 accumulates angular error. When the loop tries to close, there's a ~17° mismatch.

### Pattern 3: Positional Impossibility
```
Joint A: Part1 → Part2, AtPointConstraintIJx violation = 50.0
Joint B: Part1 → Part3, AtPointConstraintIJx violation = -50.0
```
**Diagnosis:** Part1 is being pulled 50mm in opposite directions along X. It cannot satisfy both positional constraints.

### Pattern 4: Part Too Short/Long
```
Joint A: Part1 → Part2, world_positions match
Joint B: Part1 → Part3, world_positions differ by [100, 0, 0]
```
**Diagnosis:** Part1's attachment points are 100mm apart in the model, but the assembly requires them to span a different distance.

### Pattern 5: Anti-Parallel Configuration (Now Handled)
```
MbD: Protecting 2 essential constraint(s) from removal: 15 16
```
**Note:** The solver now uses iterative constraint protection to handle anti-parallel (180°) configurations correctly. If you see "Protecting essential constraint(s)" messages followed by successful convergence, this is the solver learning which constraints are essential vs. truly redundant.

## Resolution Strategies

Based on the diagnosis, recommend one or more of:

| Strategy | When to Use |
|----------|-------------|
| **Rotate part** | Orientation violations, part at wrong angle |
| **Reposition part** | Position violations, part at wrong location |
| **Change LCS selection** | Wrong attachment point being used |
| **Resize part** | Part dimensions don't match required span |
| **Change joint type** | Over-constrained; need DOF for flexibility |
| **Remove joint** | Redundant or conflicting constraint |
| **Add intermediate part** | Gap that can't be spanned by existing parts |

## Programmatic Parsing

The diagnostic can be parsed programmatically by:
1. Looking for lines between `---BEGIN:INCONSISTENT_CONSTRAINTS---` and `---END:INCONSISTENT_CONSTRAINTS---`
2. Parsing the content as YAML

Example Python snippet:
```python
import yaml
import re

def parse_inconsistent_constraints(log_output):
    match = re.search(
        r'---BEGIN:INCONSISTENT_CONSTRAINTS---\n(.*?)---END:INCONSISTENT_CONSTRAINTS---',
        log_output,
        re.DOTALL
    )
    if match:
        yaml_content = match.group(1)
        # Skip the first line (human-readable message)
        yaml_lines = yaml_content.split('\n')[1:]
        return yaml.safe_load('\n'.join(yaml_lines))
    return None
```

## Key Principles

1. **The solver reports symptoms, not root causes.** Multiple violations may stem from a single incorrectly placed part.

2. **Start with the `affected_part`.** It appears in the most violations and is often (but not always) the misconfigured one.

3. **Use `relative_angle_degrees` for quick triage.** High values indicate orientation problems; low values suggest positional issues.

4. **Orientation errors compound.** A 10° error at one joint becomes 20° by the time a chain reaches another joint.

5. **World positions reveal ground truth.** They show where parts actually ended up in the best-effort state.

6. **The 3D view shows best-effort.** Parts are positioned at the closest achievable configuration, not their initial positions.

7. **Zero-violation constraints are filtered out.** Only constraints that actually couldn't be satisfied are reported.

8. **Iterative protection handles edge cases.** Anti-parallel configurations and configuration-dependent singularities are now handled automatically through iterative constraint protection.
