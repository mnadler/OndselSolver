# Plan: Iterative Constraint Protection for Redundancy Detection

## Problem Statement

When OndselSolver encounters a singular Jacobian matrix during Newton-Raphson iteration, it removes constraints as "redundant." However, the solver cannot distinguish between:

1. **Truly redundant constraints**: Algebraically dependent on other constraints (safe to remove)
2. **Configuration-dependent singularities**: Constraints that are essential but have a temporarily singular Jacobian (e.g., 180° anti-parallel quaternion orientations)

The current behavior removes essential constraints, causing the solver to converge to **wrong solutions** where parts are massively misaligned (e.g., 106° off instead of 0°).

## Objectives

**(a) Correct Assemblies**: Solve successfully, removing only truly redundant constraints (no information loss)

**(b) Inconsistent Assemblies**: Return the best possible visualization with descriptive error text

## Solution: Iterative Protection Approach

### Core Insight

After the solver converges, we can determine which removed constraints were essential by checking their residuals:
- **residual ≈ 0**: Constraint was truly redundant (it's satisfied even though removed)
- **residual ≠ 0**: Constraint was essential (should NOT have been removed)

### Algorithm

```
protectedEqnNos = {}  // Set of equation numbers that cannot be removed

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
                protectedEqnNos.add(essential_constraints)
                reactivate_all_constraints()
                continue  // Retry with protections

    catch SingularMatrixError(eqnNos):
        // Filter out protected constraints
        toRemove = eqnNos - protectedEqnNos

        if toRemove is empty:
            // Can't remove anything, can't converge
            throw InconsistentConstraintsError with current visual state

        remove_constraints(toRemove)
        continue
```

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
| Inconsistent assembly | Hits "can't remove, can't converge" state, reports error with visualization |

---

## Implementation Details

### File: `OndselSolver/PosICNewtonRaphson.h`

Add new member variable:

```cpp
// In class PosICNewtonRaphson
std::set<size_t> protectedEqnNos;  // Equation numbers that cannot be removed as redundant
```

### File: `OndselSolver/PosICNewtonRaphson.cpp`

#### 1. Modify `run()` - Clear protections at start

At the beginning of `run()`, clear the protection set for a fresh solve:

```cpp
void PosICNewtonRaphson::run()
{
    // Clear any previous tracking
    removedEqnNos = nullptr;
    removedRhsAtDetection = nullptr;
    protectedEqnNos.clear();

    while (true) {
        // ... existing code ...
    }
}
```

#### 2. Modify `SingularMatrixError` catch block - Filter protected constraints

Replace the existing catch block (around line 235) with filtering logic:

```cpp
catch (const SingularMatrixError& ex) {
    auto redundantEqnNos = ex.getRedundantEqnNos();
    auto rhsAtDetection = ex.getRhsValues();

    // Filter out protected constraints
    auto toRemove = std::make_shared<FullColumn<size_t>>();
    auto toRemoveRhs = std::make_shared<std::vector<double>>();

    for (size_t i = 0; i < redundantEqnNos->size(); i++) {
        size_t eqnNo = redundantEqnNos->at(i);
        if (protectedEqnNos.find(eqnNo) == protectedEqnNos.end()) {
            toRemove->push_back(eqnNo);
            if (rhsAtDetection && i < rhsAtDetection->size()) {
                toRemoveRhs->push_back(rhsAtDetection->at(i));
            }
        }
    }

    // If all singular constraints are protected, we're stuck
    if (toRemove->empty()) {
        // Reactivate constraints so diagnostic code can access them
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->reactivateRedundantConstraints();
        });

        // Build error with current state
        throw InconsistentConstraintsError(
            "Cannot solve: all singular constraints are protected",
            redundantEqnNos, rhsAtDetection);
    }

    // Remove only non-protected constraints
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        item->removeRedundantConstraints(toRemove);
    });
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        item->constraintsReport();
    });
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        item->setqsu(qsuOld);
    });

    // Store for post-convergence verification
    if (!toRemoveRhs->empty()) {
        if (!removedEqnNos) {
            removedEqnNos = std::make_shared<std::vector<size_t>>();
            removedRhsAtDetection = std::make_shared<std::vector<double>>();
        }
        for (size_t i = 0; i < toRemove->size(); i++) {
            removedEqnNos->push_back(toRemove->at(i));
            if (i < toRemoveRhs->size()) {
                removedRhsAtDetection->push_back(toRemoveRhs->at(i));
            }
        }
    }
}
```

#### 3. Modify `verifyRemovedConstraintsAtConvergence()` - Add protection logic

Update the function to protect essential constraints instead of just throwing an error:

```cpp
bool PosICNewtonRaphson::verifyRemovedConstraintsAtConvergence()
{
    // If no constraints were removed, nothing to verify
    if (!removedEqnNos || removedEqnNos->empty()) {
        return false;  // No retry needed
    }

    // Get the current constraint residuals at the converged state
    std::vector<std::pair<size_t, double>> inconsistentConstraints;
    double consistencyTolerance = 1.0e-6;

    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        auto joint = std::dynamic_pointer_cast<Joint>(item);
        if (joint) {
            joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
                if (con->isRedundant()) {
                    auto redunCon = std::static_pointer_cast<RedundantConstraint>(con);
                    auto wrappedCon = redunCon->constraint;
                    wrappedCon->postPosICIteration();
                    double residual = wrappedCon->aG;
                    if (std::abs(residual) > consistencyTolerance) {
                        inconsistentConstraints.push_back({wrappedCon->iG, residual});
                    }
                }
            });
        }
    });

    // Also check Part constraints (aGeu, aGabs)
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        auto part = std::dynamic_pointer_cast<Part>(item);
        if (part && part->partFrame) {
            auto aGeu = part->partFrame->aGeu;
            if (aGeu && aGeu->isRedundant()) {
                auto redunCon = std::static_pointer_cast<RedundantConstraint>(aGeu);
                auto wrappedCon = redunCon->constraint;
                wrappedCon->postPosICIteration();
                double residual = wrappedCon->aG;
                if (std::abs(residual) > consistencyTolerance) {
                    inconsistentConstraints.push_back({wrappedCon->iG, residual});
                }
            }
            if (part->partFrame->aGabs) {
                for (auto& aGab : *(part->partFrame->aGabs)) {
                    if (aGab->isRedundant()) {
                        auto redunCon = std::static_pointer_cast<RedundantConstraint>(aGab);
                        auto wrappedCon = redunCon->constraint;
                        wrappedCon->postPosICIteration();
                        double residual = wrappedCon->aG;
                        if (std::abs(residual) > consistencyTolerance) {
                            inconsistentConstraints.push_back({wrappedCon->iG, residual});
                        }
                    }
                }
            }
        }
    });

    if (!inconsistentConstraints.empty()) {
        // Protect these constraints from future removal
        std::ostringstream protectedMsg;
        protectedMsg << "MbD: Protecting " << inconsistentConstraints.size()
                     << " essential constraint(s) from removal:";
        for (const auto& pair : inconsistentConstraints) {
            protectedEqnNos.insert(pair.first);
            protectedMsg << " " << pair.first;
        }
        system->logString(protectedMsg.str());

        // Reactivate ALL constraints and retry
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->reactivateRedundantConstraints();
        });

        // CRITICAL: Reset positions to initial state before retry
        // Without this, we start from the wrong converged state and hit different singularities
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->setqsu(qsuOld);
        });

        // Clear tracking for fresh retry
        removedEqnNos = nullptr;
        removedRhsAtDetection = nullptr;

        return true;  // Signal to retry
    }

    return false;  // No retry needed - all removed constraints are truly redundant
}
```

### Required Include

Ensure `<set>` is included in `PosICNewtonRaphson.cpp`:

```cpp
#include <set>
```

---

## Code Flow After Changes

```
run()
├── protectedEqnNos.clear()
├── while (true)
│   ├── try
│   │   ├── preRun() → "Assembling system"
│   │   ├── iterate() → Newton-Raphson iterations
│   │   │   └── solveEquations() may throw SingularMatrixError
│   │   ├── verifyRemovedConstraintsAtConvergence()
│   │   │   ├── If essential constraints found:
│   │   │   │   ├── Add to protectedEqnNos
│   │   │   │   ├── reactivateRedundantConstraints()
│   │   │   │   ├── setqsu(qsuOld) ← CRITICAL: reset to initial positions
│   │   │   │   └── return true → continue loop
│   │   │   └── If all removed are truly redundant:
│   │   │       └── return false → break loop (SUCCESS)
│   │   └── break (SUCCESS)
│   │
│   └── catch SingularMatrixError
│       ├── Filter out protectedEqnNos from redundantEqnNos
│       ├── If toRemove is empty:
│       │   └── throw InconsistentConstraintsError (STUCK)
│       ├── removeRedundantConstraints(toRemove)
│       └── continue loop
```

---

## Test Cases

### 1. Flat/Exploded Assembly (Correct Constraints)
- **Input**: Assembly with all parts at origin, anti-parallel orientations (180°)
- **Expected**:
  - First iterations: removes some quaternion constraints
  - At convergence: detects these were essential (large residuals)
  - Protects them, reactivates all, retries
  - Eventually converges correctly with only truly redundant constraints removed
- **Verification**: All joints have `relative_angle_degrees ≈ 0`

### 2. Assembly with True Redundancy
- **Input**: Over-constrained mechanism with algebraically redundant constraints
- **Expected**: Converges on first try, redundant constraints have residual ≈ 0
- **Verification**: Removed constraints stay removed, assembly correct

### 3. Inconsistent Assembly
- **Input**: Assembly with contradictory constraints (e.g., x=0 AND x=1)
- **Expected**:
  - Iteratively protects constraints
  - Eventually all singular constraints are protected
  - Throws `InconsistentConstraintsError` with diagnostic
- **Verification**: Error message contains joint diagnostics

### 4. Regression Tests
- Run existing test suite to ensure no regressions
- Specifically test:
  - `correct assembly perturbation` tests
  - Any tests with redundant constraints

---

## Complexity Analysis

- **Worst case**: O(N) iterations where N = number of constraints
- **Typical case**: 1-3 iterations for most assemblies
- **Each iteration**: Full Newton-Raphson solve

---

## Removed Code

The following code from the previous perturbation-based implementation attempt was removed:

1. **`hasRetriedWithPerturbation`** - No longer using perturbation approach
2. **Perturbation logic** - Quaternion perturbation code that tried to escape singularities
3. **`quaternionAntiParallelThreshold`** - No longer classifying by threshold

The iterative protection approach supersedes the perturbation approach entirely.

---

## Files to Modify

| File | Changes |
|------|---------|
| `OndselSolver/PosICNewtonRaphson.h` | Add `std::set<size_t> protectedEqnNos` member |
| `OndselSolver/PosICNewtonRaphson.cpp` | Modify `run()`, `SingularMatrixError` catch, `verifyRemovedConstraintsAtConvergence()` |

---

## Summary

This approach:
1. **Uses existing infrastructure** (`removeRedundantConstraints`, `reactivateRedundantConstraints`)
2. **Learns empirically** which constraints are essential (no a priori classification needed)
3. **Guarantees termination** (finite constraints → finite iterations)
4. **Handles all cases**: correct assemblies, anti-parallel configurations, and inconsistent assemblies
5. **Minimal code changes** (only PosICNewtonRaphson modified)
