# 180° Singularity Correction: Analysis and Fix

## Executive Summary

The OndselSolver's assembly solver was experiencing infinite loops when solving complex assemblies. The root cause was the 180° singularity correction mechanism fighting against itself: each correction escaped one singularity but created another, and the solver kept "correcting" the same parts in an endless cycle.

**The Fix:** Correct each part at most once.

**Why This Is Correct:** A 180° singularity correction is fundamentally a one-time escape hatch. If the true solution doesn't require 180° orientation, one nudge is sufficient. If it does require 180°, no amount of nudging will help—the assembly is geometrically singular and should fail cleanly rather than loop forever.

---

## Problem Statement

When solving assemblies with parts in near-180° orientations, the solver would timeout instead of converging. Debug logging revealed the solver was oscillating: correcting Part A about the X-axis, then about the Z-axis, then X again, then Z again, indefinitely.

**Observed behavior in logs:**
```
19:11:48 setSize=4   ← Corrections accumulated
19:11:51 setSize=0   ← New solver instance, memory lost
19:11:57 setSize=0   ← Another new instance
19:12:04 setSize=0   ← And again... until timeout
```

The solver was losing memory of which parts it had already corrected, causing it to re-correct the same parts repeatedly.

---

## Context

### What the 180° Singularity Is

When two coordinate frames are oriented exactly 180° apart (anti-parallel), the rotation axis between them becomes mathematically undefined. The Jacobian matrix becomes singular, and Newton-Raphson iteration cannot proceed.

### What the Correction Does

The solver detects this condition and applies a small rotation (a few degrees) to nudge the part away from the degenerate 180° configuration. This allows Newton-Raphson to continue. The solver tracks corrected parts in a set called `correctedPartsFor180` to avoid re-correcting them.

### Why the Memory Was Being Lost

The solver has a retry mechanism. After each solve attempt, it verifies that constraints marked as "redundant" are truly redundant. If not, it reactivates them and retries:

```cpp
runPosIC();                      // First attempt
while (needToRedoPosIC()) {      // Verify constraints
    runPosIC();                  // Retry if needed
}
```

The problem was in `runPosIC()`:

```cpp
void SystemSolver::runPosIC() {
    icTypeSolver = CREATE<PosICNewtonRaphson>::With();  // NEW instance every time
    icTypeSolver->run();
}
```

Every retry created a new solver instance. Since `correctedPartsFor180` is a member variable, each new instance started with an empty set—no memory of previous corrections.

### The Oscillation Mechanism

1. **Attempt 1:** Part A at 180° on X-axis → correct about X → `set = {A}`
2. **Retry triggered:** Constraint verification requires another attempt
3. **Attempt 2:** New solver, `set = {}` → Part A now at 180° on Z-axis (caused by X correction) → correct about Z
4. **Retry triggered**
5. **Attempt 3:** New solver, `set = {}` → Part A back at 180° on X-axis (caused by Z correction) → correct about X
6. **Forever:** X → Z → X → Z → timeout

---

## The Fix

Two changes ensure correction memory persists:

### Change 1: Reuse Solver Instance

**File:** `OndselSolver/SystemSolver.cpp`

```cpp
void SystemSolver::runPosIC()
{
    // Only create new solver if we don't have a PosICNewtonRaphson already.
    // Reusing preserves correctedPartsFor180 across retries, preventing
    // the same part from being 180°-corrected multiple times (oscillation).
    if (!std::dynamic_pointer_cast<PosICNewtonRaphson>(icTypeSolver)) {
        icTypeSolver = CREATE<PosICNewtonRaphson>::With();
    }
    icTypeSolver->setSystem(this);
    icTypeSolver->run();
}
```

### Change 2: Don't Clear the Set

**File:** `OndselSolver/PosICNewtonRaphson.cpp`

Remove `correctedPartsFor180.clear()` from the `run()` method. The set now accumulates corrections across the entire solve operation.

```cpp
void PosICNewtonRaphson::run()
{
    // Clear per-run state
    removedEqnNos = nullptr;
    removedRhsAtDetection = nullptr;
    protectedConstraints.clear();
    allRemovedConstraints.clear();
    // NOTE: correctedPartsFor180 is intentionally NOT cleared here.
    // See docs/fix-180-degree-correction-oscillation.md for rationale.

    // ... rest of run() ...
}
```

---

## Discussion: Why This Fix Is Correct Across All Cases

### The Fundamental Insight

**A 180° correction is a one-time escape hatch, not an iterative solution strategy.**

There are only two possible scenarios:

| Scenario | What Happens | Outcome |
|----------|--------------|---------|
| True solution is NOT at 180° | One correction nudges part away from singularity, solver converges to the true solution | Success |
| True solution IS at 180° | Correction nudges away, solver tries to return to 180°, hits singularity again | Should fail cleanly |

In the second scenario, re-correcting doesn't help—it just delays failure while burning cycles. The correction and the constraints are fundamentally incompatible. The assembly requires a 180° relationship that is mathematically singular.

**The fix makes the solver fail fast and cleanly rather than loop forever.**

### What About Multiple Parts?

The fix only prevents re-correction of the *same* part. Different parts can still be corrected independently:

- Part A hits singularity → corrected → `set = {A}`
- Part B hits singularity → not in set → corrected → `set = {A, B}`

This is correct behavior.

### What About BFS Incremental Placement?

FreeCAD uses BFS incremental placement, calling `runPreDrag()` for each part being placed. Each call creates a new `System` object, which resets the correction set. This is correct—each BFS step places different parts with independent geometry. The oscillation problem occurs *within* a single step's Newton-Raphson iteration, not between steps.

### Potential Counterexample (And Why It's Invalid)

One might ask: "What if constraint reactivation legitimately requires a different correction axis?"

If Part A needs correction on a different axis after reactivation, it means the true solution (with all constraints active) requires Part A at 180° on that axis. This is the same fundamental problem—we'd just oscillate on a different axis. Re-correction wouldn't lead to convergence; it would lead to a different oscillation pattern.

**The counterexample is actually another case where the assembly is geometrically problematic and should fail.**

---

## Safety Analysis

### Between Different `solve()` Calls

Each `solve()` call creates a completely new object hierarchy:

```
assembly.solve()
  → ASMTAssembly::runKINEMATIC()
    → mbdSystem = std::make_shared<System>()     ← NEW System
      → systemSolver = std::make_shared<SystemSolver>()  ← NEW SystemSolver
        → icTypeSolver = nullptr  ← Fresh state
```

Fresh solver with empty `correctedPartsFor180` for each new solve operation.

### Within a `solve()` Call

We reuse the solver across retries, preserving `correctedPartsFor180`. Corrections persist until the solve operation completes.

### After Velocity/Acceleration Solving

The code calls `runVelIC()` and `runAccIC()` after position solving, which overwrite `icTypeSolver` with different solver types. If there's an outer loop iteration, `dynamic_pointer_cast<PosICNewtonRaphson>` returns null, creating a fresh position solver.

---

## Files Changed

1. `OndselSolver/SystemSolver.cpp` - Reuse `PosICNewtonRaphson` instance in `runPosIC()`
2. `OndselSolver/PosICNewtonRaphson.cpp` - Remove `correctedPartsFor180.clear()` from `run()`

---

## Testing

To verify the fix:

1. Build the solver with debug logging enabled
2. Run a complex assembly that previously caused timeout
3. Confirm in logs that:
   - `setSize` increases across retries (not resetting to 0)
   - "SKIPPING - part already in set" messages appear for repeated corrections
   - Assembly solves successfully without timeout

---

## Future Considerations

### Enhanced Error Reporting

When the solver fails to converge after corrections have been applied, consider reporting a clear error message:

> "Assembly contains 180° singularity that cannot be resolved. Check that connected parts are not required to be exactly anti-parallel."

This would help users understand that the problem is with the assembly design, not the solver.
