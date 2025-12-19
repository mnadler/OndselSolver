# Iterative Constraint Protection for Redundancy Detection

## Problem Statement

When OndselSolver encounters a singular Jacobian matrix during Newton-Raphson iteration, it removes constraints as "redundant." However, the solver cannot distinguish between:

1. **Truly redundant constraints**: Algebraically dependent on other constraints (safe to remove)
2. **Configuration-dependent singularities**: Constraints that are essential but have a temporarily singular Jacobian (e.g., 180° anti-parallel quaternion orientations)

The previous behavior removed essential constraints, causing the solver to converge to **wrong solutions** where parts are massively misaligned (e.g., 106° off instead of 0°).

## Objectives

**(a) Correct Assemblies**: Solve successfully, removing only truly redundant constraints (no information loss)

**(b) Inconsistent Assemblies**: Display the best possible visualization with descriptive error text in the console

## Solution: Iterative Protection Approach

### Core Insight

After the solver converges, we can determine which removed constraints were essential by checking their residuals:
- **residual ≈ 0**: Constraint was truly redundant (it's satisfied even though removed)
- **residual ≠ 0**: Constraint was essential (should NOT have been removed)

### Algorithm

```
protectedConstraints = {}  // Set of Constraint* that cannot be removed
bestEffortState = nullptr  // Saved converged state before retry

while true:
    try:
        run Newton-Raphson solver

        if converged:
            essential_constraints = check_removed_constraints_at_convergence()

            if essential_constraints is empty:
                SUCCESS - all removed constraints were truly redundant
                break
            else:
                // Some removed constraints were essential
                bestEffortState = save_current_positions()
                protectedConstraints.add(essential_constraints)
                reactivate_all_constraints()
                reset_to_initial_positions()
                continue  // Retry with protections

    catch SingularMatrixError(eqnNos):
        // Filter out protected constraints
        toRemove = eqnNos - protectedConstraints

        if toRemove is empty:
            // Can't remove anything, can't converge
            restore bestEffortState (if available)
            log YAML diagnostic
            return normally  // Don't throw - let FreeCAD display best-effort state

        remove_constraints(toRemove)
        continue
```

### Key Design Decisions

1. **Using Constraint pointers instead of equation numbers**: Equation numbers (`iG`) change between retries when constraints are reactivated. Using `std::set<Constraint*>` ensures protected constraints remain protected across retries.

2. **Saving best-effort state**: Before retrying with protections, we save the current converged state. If the assembly is ultimately inconsistent, this best-effort state is restored for visualization.

3. **Returning instead of throwing**: When constraints are inconsistent, instead of throwing `InconsistentConstraintsError`, we log the diagnostic and return normally. This allows FreeCAD to call `setNewPlacements()` and display the best-effort state.

### Termination Guarantee

Each iteration either:
1. **Succeeds** (no essential constraints incorrectly removed), OR
2. **Protects at least one more constraint**

Since there are finite constraints, the algorithm terminates in at most N iterations (where N = number of constraints).

### Outcomes

| Scenario | Result |
|----------|--------|
| Correct assembly with redundancy | Converges correctly, only truly redundant constraints removed |
| Correct assembly with anti-parallel parts | Iteratively learns which constraints are essential, eventually converges |
| Inconsistent assembly | Returns best-effort state with YAML diagnostic logged to console |

---

## Implementation

### File: `OndselSolver/PosICNewtonRaphson.h`

```cpp
class PosICNewtonRaphson : public AnyPosICNewtonRaphson
{
public:
    void run() override;
    void preRun() override;
    void assignEquationNumbers() override;
    bool isConverged() override;
    void handleSingularMatrix() override;
    void lookForRedundantConstraints();
    bool verifyRemovedConstraintsAtConvergence();  // Returns true if retry needed

    std::shared_ptr<std::vector<size_t>> pivotRowLimits;

    // Track constraints removed as potentially-redundant for post-convergence verification
    std::shared_ptr<std::vector<size_t>> removedEqnNos;
    std::shared_ptr<std::vector<double>> removedRhsAtDetection;

    // Constraint pointers that cannot be removed as redundant (learned to be essential)
    // Using pointers instead of equation numbers because iG changes on each retry
    std::set<Constraint*> protectedConstraints;

    // Track ALL constraints ever removed across all iterations for final violation reporting
    std::set<Constraint*> allRemovedConstraints;

    // Best effort converged state - saved before retrying, restored before returning on error
    FColDsptr bestEffortState;
};
```

### File: `OndselSolver/PosICNewtonRaphson.cpp`

#### Structure of `run()`

```cpp
void PosICNewtonRaphson::run()
{
    // Clear any previous tracking
    removedEqnNos = nullptr;
    removedRhsAtDetection = nullptr;
    protectedConstraints.clear();
    allRemovedConstraints.clear();

    try {  // OUTER try - catches ALL InconsistentConstraintsError for YAML processing
        while (true) {
            try {
                preRun();
                initializeLocally();
                initializeGlobally();
                iterate();
                postRun();

                // After successful convergence, verify removed constraints are satisfied
                if (verifyRemovedConstraintsAtConvergence()) {
                    continue;  // Retry the solve with protected constraints
                }
                break;  // Success
            }
            catch (const SingularMatrixError& ex) {
                // Filter out protected constraints from removal candidates
                // If all are protected, throw InconsistentConstraintsError
                // Otherwise, remove non-protected constraints and retry
                // ...
            }
        }
    }
    catch (InconsistentConstraintsError& ex) {
        // Set bestEffortState (saved before retry in verifyRemovedConstraintsAtConvergence)
        // Update Parts with postPosICIteration() to refresh marker frame positions
        // Build YAML diagnostic with joint/constraint details at best-effort state
        // Log diagnostic and call updateFromMbD() to push to FreeCAD
        system->logString(diagnostic);
        return;  // Let FreeCAD display best-effort state
    }
}
```

#### `verifyRemovedConstraintsAtConvergence()`

This function checks if any constraints marked as redundant actually have non-zero residuals at convergence:

```cpp
bool PosICNewtonRaphson::verifyRemovedConstraintsAtConvergence()
{
    // If no constraints were removed, nothing to verify
    if (!removedEqnNos || removedEqnNos->empty()) {
        return false;
    }

    // Check residuals of all RedundantConstraint wrappers
    std::map<Constraint*, double> essentialConstraints;
    double consistencyTolerance = 1.0e-6;

    // Check joint constraints...
    // Check part constraints (aGeu, aGabs)...

    if (!essentialConstraints.empty()) {
        // Log which constraints are being protected
        // Save bestEffortState before retry
        bestEffortState = std::make_shared<FullColumn<double>>(nqsu);
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->fillqsu(bestEffortState);
        });

        // Add to protected set
        for (const auto& pair : essentialConstraints) {
            protectedConstraints.insert(pair.first);
        }

        // Reactivate ALL constraints and reset to initial positions
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->reactivateRedundantConstraints();
        });
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->setqsu(qsuOld);
        });

        // Clear tracking for fresh retry
        removedEqnNos = nullptr;
        removedRhsAtDetection = nullptr;

        return true;  // Signal to retry
    }

    return false;  // All removed constraints are truly redundant
}
```

### File: `OndselSolver/ASMTAssembly.cpp`

Since we no longer throw exceptions for inconsistent constraints, `runKINEMATIC()` is simplified:

```cpp
void MbD::ASMTAssembly::runKINEMATIC()
{
    mbdSystem = std::make_shared<System>();
    mbdSystem->externalSystem->asmtAssembly = this;
    mbdSystem->runKINEMATIC(mbdSystem);
}
```

---

## YAML Diagnostic Format

When constraints are inconsistent, a YAML diagnostic is logged to the console:

```yaml
---BEGIN:INCONSISTENT_CONSTRAINTS---
Constraints are geometrically inconsistent (no solution exists)
  affected_part: "PartName"
  total_violation: 0.123456
  inconsistent_equation_count: 3
  nqsu: 14
  joints:
    - name: "JointName"
      type: "RevoluteJoint"
      part_i: "Part1"
      part_j: "Part2"
      lcs_i:
        name: "LCS_Origin"
        position_on_part: [0.0, 0.0, 0.0]
        world_position: [1.0, 2.0, 3.0]
        world_quaternion: [0.0, 0.0, 0.0, 1.0]
      lcs_j:
        name: "LCS_Origin"
        position_on_part: [0.0, 0.0, 0.0]
        world_position: [1.1, 2.0, 3.0]
        world_quaternion: [0.0, 0.0, 0.0, 1.0]
      relative_angle_degrees: 0.5
      relative_quaternion: [0.0, 0.0, 0.004, 1.0]
      inconsistent_constraints:
        - type: "TranslationConstraintIJ directioni"
          violation: 0.1
  part_constraints:
    - part: "Part1"
      type: "EulerParameterConstraint"
      equation_number: 7
      violation: 0.0001
---END:INCONSISTENT_CONSTRAINTS---
```

The `---BEGIN:` and `---END:` markers allow Python code to parse the diagnostic if needed.

### Ground Joint Filtering

Ground joints (joints that fix a part to the assembly/world) are automatically filtered from the diagnostic output because:

1. **Ground joints define the reference frame** - they cannot be "inconsistent" by themselves
2. **If a ground joint appears violated**, it means OTHER joints are pulling the part away from where it was grounded - the conflict is with those other joints, not the ground joint

**Detection**: Ground parts are identified by their path structure:
- Ground/assembly: `/OndselAssembly` (no nested path)
- Regular part: `/OndselAssembly/Part#Link` (has nested path with second `/`)

The diagnostic only shows joints between actual parts, not joints connecting parts to the assembly origin.

---

## Code Flow

```
run()
├── Clear tracking (protectedConstraints, allRemovedConstraints, etc.)
├── OUTER try
│   ├── while (true)
│   │   ├── INNER try
│   │   │   ├── preRun() → "Assembling system"
│   │   │   ├── iterate() → Newton-Raphson iterations
│   │   │   │   └── solveEquations() may throw SingularMatrixError
│   │   │   ├── verifyRemovedConstraintsAtConvergence()
│   │   │   │   ├── If essential constraints found:
│   │   │   │   │   ├── Save bestEffortState
│   │   │   │   │   ├── Add to protectedConstraints
│   │   │   │   │   ├── reactivateRedundantConstraints()
│   │   │   │   │   ├── setqsu(qsuOld) ← reset to initial positions
│   │   │   │   │   └── return true → continue loop
│   │   │   │   └── If all removed are truly redundant:
│   │   │   │       └── return false → break loop (SUCCESS)
│   │   │   └── break (SUCCESS)
│   │   │
│   │   └── catch SingularMatrixError
│   │       ├── Filter out protectedConstraints
│   │       ├── If toRemove is empty:
│   │       │   └── throw InconsistentConstraintsError
│   │       ├── removeRedundantConstraints(toRemove)
│   │       ├── Track in allRemovedConstraints
│   │       └── continue loop
│   │
│   └── catch InconsistentConstraintsError
│       ├── Reactivate constraints
│       ├── Set bestEffortState (saved in verifyRemovedConstraintsAtConvergence before retry)
│       ├── Update Parts with postPosICIteration() ← updates marker frame world positions
│       ├── Compute constraint residuals at best-effort state
│       ├── Build YAML diagnostic (uses best-effort positions)
│       ├── Log diagnostic
│       ├── Call updateFromMbD() ← pushes best-effort state to FreeCAD
│       └── RETURN (don't throw) ← FreeCAD displays best-effort state
```

**Note:** Both the YAML diagnostic and the FreeCAD display use `bestEffortState` - the solver's best converged configuration before determining inconsistency. This allows users to see where the solver got stuck, which is more useful for debugging than showing initial positions.

---

## Test Cases

### 1. Flat/Exploded Assembly (Correct Constraints)
- **Input**: Assembly with all parts at origin, anti-parallel orientations (180°)
- **Expected**:
  - First iterations: removes some quaternion constraints
  - At convergence: detects these were essential (large residuals)
  - Saves bestEffortState, protects them, reactivates all, retries
  - Eventually converges correctly with only truly redundant constraints removed
- **Verification**: All joints have `relative_angle_degrees ≈ 0`

### 2. Assembly with True Redundancy
- **Input**: Over-constrained mechanism with algebraically redundant constraints
- **Expected**: Converges on first try, redundant constraints have residual ≈ 0
- **Verification**: Removed constraints stay removed, assembly correct

### 3. Inconsistent Assembly
- **Input**: Assembly with contradictory constraints (e.g., parts that can't physically connect)
- **Expected**:
  - Iteratively protects constraints
  - Eventually all singular constraints are protected
  - Logs YAML diagnostic and returns normally
  - FreeCAD displays best-effort state
- **Verification**: Console shows YAML diagnostic, 3D view shows best-effort positions

---

## Files Modified

| File | Changes |
|------|---------|
| `OndselSolver/PosICNewtonRaphson.h` | Added `protectedConstraints`, `allRemovedConstraints`, `bestEffortState` members |
| `OndselSolver/PosICNewtonRaphson.cpp` | Implemented iterative protection in `run()`, `verifyRemovedConstraintsAtConvergence()`, filtering in `SingularMatrixError` catch, ground joint filtering in YAML diagnostic output |
| `OndselSolver/ASMTAssembly.cpp` | Simplified `runKINEMATIC()` (removed try-catch since we don't throw) |

---

## Summary

This implementation:
1. **Uses existing infrastructure** (`removeRedundantConstraints`, `reactivateRedundantConstraints`)
2. **Learns empirically** which constraints are essential (no a priori classification needed)
3. **Guarantees termination** (finite constraints → finite iterations)
4. **Handles all cases**: correct assemblies, anti-parallel configurations, and inconsistent assemblies
5. **Displays best-effort state** for inconsistent assemblies (returns normally instead of throwing)
6. **Provides detailed diagnostics** via YAML format in console log
7. **Filters ground joints** from diagnostics (they define the reference frame and can't be truly inconsistent)
