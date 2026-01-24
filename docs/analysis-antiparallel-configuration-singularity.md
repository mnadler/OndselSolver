# Analysis: Anti-Parallel (180°) Configuration Singularity in Quaternion Constraints

## Executive Summary

This document analyzes the failure mode observed when the OndselSolver attempts to solve a simple two-part assembly where the Local Coordinate Systems (LCS) are oriented 180° apart (anti-parallel). The solver incorrectly removes quaternion constraints as "redundant" when they are in fact essential, leading to a failure to converge.

**Key Finding**: This is a **configuration singularity**, not a constraint redundancy. The Jacobian matrix becomes rank-deficient at exactly 180° rotation, but the constraints themselves remain essential. The current solver's redundancy detection mechanism conflates these two distinct mathematical conditions.

---

## 1. Problem Description

### 1.1 Observed Behavior

When two parts have their LCS frames oriented 180° apart:
1. The solver detects a singular Jacobian matrix during Newton-Raphson iteration
2. It incorrectly interprets this as a redundant constraint (e.g., `QuaternionConstraintx`)
3. The constraint is removed and wrapped as `RedundantConstraint`
4. After convergence, `verifyRemovedConstraintsAtConvergence()` detects a large violation (1.0)
5. The constraint is marked as "protected" and the solve is retried
6. On retry, positions are reset to `qsuOld` (the original 180° configuration)
7. The matrix is still singular at 180°, but the constraint is now protected
8. Error: "all singular constraints are protected"

### 1.2 Solver Log Evidence

```
MbD: Checking for redundant constraints.
GESpMatFullPvPosIC: Pivot too small at row 9, treating as redundant
Removing constraint equation 9 (QuaternionConstraintx) as redundant
...
MbD: Protecting 1 essential constraint(s) from removal: 9
...
Cannot solve: all singular constraints are protected (geometrically inconsistent)
  relative_angle_degrees: 180.000000
  violation: 1.000000
```

---

## 2. Mathematical Analysis

### 2.1 Quaternion Constraint Formulation

The OndselSolver uses quaternion constraints to enforce relative orientation. For a fixed joint, three constraints enforce that the imaginary part of the relative quaternion is zero:

```
C_k = Im_k(conj(q_I) * q_J) = 0,  k ∈ {x, y, z}
```

Where:
- `q_I` = world quaternion of frame I (accounts for marker frame)
- `q_J` = world quaternion of frame J (accounts for marker frame)
- `conj(q)` = quaternion conjugate: `(-q_x, -q_y, -q_z, q_w)`

The constraint value is computed in `QuaternionIecJec.cpp` using the Hamilton product formula:
```cpp
// p = conj(qI) * qJ
// p[0] = qI3*qJ0 - qI0*qJ3 - qI1*qJ2 + qI2*qJ1  (x imaginary)
// p[1] = qI3*qJ1 + qI0*qJ2 - qI1*qJ3 - qI2*qJ0  (y imaginary)
// p[2] = qI3*qJ2 - qI0*qJ1 + qI1*qJ0 - qI2*qJ3  (z imaginary)
```

### 2.2 Why Quaternions Have a Unique Solution

Unlike perpendicular constraints (which have ambiguity at ±90°), quaternion constraints have exactly one correct orientation:
- The target is `conj(q_I) * q_J = (0, 0, 0, 1)` (identity quaternion)
- This corresponds to zero relative rotation between frames
- There is no second solution because quaternions encode orientation directly, not just direction cosines

**Reference**: [RBDL Joint Modeling](https://rbdl.github.io/df/dbe/joint_description.html): "When using quaternion joint type the model does not suffer from singularities [in the representation sense]. A quaternion has four values but represents only three degrees of freedom."

### 2.3 The 180° Configuration Singularity

At exactly 180° relative rotation, the quaternion takes the form:
- `q_rel = conj(q_I) * q_J = (±1, 0, 0, 0)` or permutations

The constraint values are:
- `C_x = ±1` (or 0)
- `C_y = 0` (or ±1)
- `C_z = 0` (or ±1)

The issue is in the **Jacobian**, not the constraint values. At 180°:

```
∂C_k/∂q = function of (q_I, q_J)
```

When `q_I` and `q_J` differ by 180°, the Jacobians of all three constraints become **colinear** or **zero** in certain components. This is because all three constraints are "trying to say the same thing": rotate by 180° in some direction.

**Key Insight**: The Jacobian rank deficiency does NOT mean the constraints are redundant. It means the linearized system cannot distinguish which direction to rotate to fix the 180° error. This is a **configuration singularity**, distinct from **constraint redundancy**.

### 2.4 Configuration Singularity vs. Constraint Redundancy

| Property | Configuration Singularity | Constraint Redundancy |
|----------|--------------------------|----------------------|
| Cause | Jacobian rank-deficient at specific pose | Multiple constraints enforce same DOF |
| Solution | Move to non-singular configuration | Remove truly redundant constraints |
| Constraint value at issue | Non-zero (constraint violated) | Zero (constraint satisfied by others) |
| Permanent? | No - singularity goes away at other configurations | Yes - constraint is always redundant |

**Reference**: [ANSYS Motion Theory](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/motion_theory/motion_theory_fundamentals_constraints.html): "Singular configuration in a multibody system occurs when the constraint equations lead to a Jacobian matrix that is linearly dependent. These conditions present configurations where the system's normal operation is disrupted... This occurs when the direction of a constraint coincides with the direction of the lost degree of freedom."

---

## 3. Literature Review: Handling Singular Jacobians

### 3.1 Levenberg-Marquardt Algorithm

The Levenberg-Marquardt (LM) algorithm is the standard solution for Newton-Raphson with singular or ill-conditioned Jacobians:

**Standard Newton-Raphson**:
```
J^T J · δ = -J^T r
```

**Levenberg-Marquardt**:
```
(J^T J + λI) · δ = -J^T r
```

Where λ is a damping parameter that:
- Regularizes the matrix, ensuring it is always invertible
- When λ → 0: Gauss-Newton behavior (fast convergence near solution)
- When λ → ∞: Gradient descent behavior (stable from poor starting points)

**Reference**: [Levenberg-Marquardt Wikipedia](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm): "The solution was suggested by Levenberg and Marquardt using an iterative sequence where a positive regularization or damping parameter is specifically added to remedy the possibility of the Jacobian matrix becoming rank-deficient."

**Reference**: [Duke LM Tutorial](https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf): "The LMA is more robust than the GNA, which means that in many cases it finds a solution even if it starts very far off the final minimum."

### 3.2 Trust-Region Methods

Trust-region methods limit the step size, which naturally handles singularities:
- Compute the Newton direction if possible
- Limit the step to stay within a "trust region" where the linearization is valid
- At singular configurations, the trust region prevents unbounded steps

### 3.3 Continuation/Homotopy Methods

For assemblies that start far from the solution:
- Define a path from the current configuration to the target
- Solve incrementally along the path
- Each intermediate step is a small perturbation, avoiding singular jumps

### 3.4 Perturbation from Singular Configurations

When detecting a singular configuration:
- Apply a small random perturbation to the current state
- The perturbation moves the system away from the singularity
- Continue iteration from the perturbed state

**Caution**: This approach is noted in literature as potentially "fragile" for general use because:
- The perturbation direction may be wrong, leading to longer convergence
- In some geometries, perturbation may push toward a different solution branch
- It requires careful magnitude selection

---

## 4. Analysis of Current Solver Behavior

### 4.1 The Core Bug

The solver's `GESpMatFullPvPosIC` performs LU decomposition with partial pivoting. When a pivot is too small (below tolerance), it reports the row as "redundant":

```cpp
// In lookForRedundantConstraints():
posICsolver->solvewithsaveOriginal(pypx, y->negated(), false);
// Throws SingularMatrixError with redundantEqnNos
```

The catch handler in `PosICNewtonRaphson::run()` then:
1. Removes the constraints with those equation numbers
2. Wraps them in `RedundantConstraint`
3. Continues iteration

**The bug**: A small pivot does NOT distinguish between:
- **True redundancy**: Constraint is linearly dependent on others (violation = 0 at solution)
- **Configuration singularity**: Jacobian is rank-deficient at current pose (violation ≠ 0)

### 4.2 Why Protected Constraints Logic Partially Works

The `verifyRemovedConstraintsAtConvergence()` function attempts to detect this case:
1. After convergence, it checks each removed constraint's actual residual
2. If |residual| > tolerance (1e-6), the constraint was essential, not redundant
3. It protects that constraint and triggers a retry

This is correct logic for **detecting** the error post-hoc. However, the retry mechanism fails because:
1. It resets to `qsuOld` (the original 180° configuration)
2. At 180°, the matrix is still singular
3. The protected constraint can't be removed
4. Result: "all singular constraints are protected"

### 4.3 Why Resetting to qsuOld Fails

The line in `verifyRemovedConstraintsAtConvergence()`:
```cpp
system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
    item->setqsu(qsuOld);  // Reset to initial state
});
```

This resets positions to the **initial configuration** before retrying. If the initial configuration was at 180°, we're back in the singular configuration. The matrix will be singular again, but now we can't remove the protected constraint.

---

## 5. Proposed Solutions

### 5.1 Option A: Levenberg-Marquardt Damping (Recommended)

**Rationale**: LM is the industry-standard solution for singular Jacobians. It addresses the root cause by making the matrix always invertible.

**Implementation**:
```cpp
// In solveEquations():
// Instead of: J^T J · δ = -J^T r
// Use: (J^T J + λI) · δ = -J^T r

// With adaptive λ:
// - Start with small λ (e.g., 1e-6)
// - If step increases error, increase λ (toward gradient descent)
// - If step decreases error, decrease λ (toward Gauss-Newton)
```

**Advantages**:
- Mathematically principled
- Well-understood convergence properties
- No special-casing for specific singularities
- Works for all configuration singularities, not just 180°

**Disadvantages**:
- More significant code change
- Requires tuning of λ adaptation strategy

### 5.2 Option B: Distinguish Singularity from Redundancy at Detection

**Rationale**: Check constraint residual at the moment singularity is detected, not just at convergence.

**Implementation**:
```cpp
// In SingularMatrixError handler:
for (size_t eqnNo : redundantEqnNos) {
    Constraint* con = findConstraintByEqnNo(eqnNo);
    double residual = con->aG;  // Current violation
    if (std::abs(residual) > 1e-4) {
        // Large violation = configuration singularity, NOT redundancy
        // Do NOT remove this constraint
        // Instead, apply perturbation or switch to LM
    }
}
```

**Advantages**:
- Fixes the immediate bug
- Relatively small code change
- Preserves existing redundancy detection for true redundancies

**Disadvantages**:
- Requires additional handling for what to do after detecting singularity
- Still needs fallback mechanism (perturbation or LM)

### 5.3 Option C: Perturbation on Retry (Minimal Change)

**Rationale**: When protected constraints prevent removal, perturb away from singularity instead of staying at qsuOld.

**Implementation**:
```cpp
// In verifyRemovedConstraintsAtConvergence(), before retry:
if (!essentialConstraints.empty()) {
    // Instead of: item->setqsu(qsuOld);
    // Apply small random perturbation to qsuOld
    auto qsuPerturbed = perturb(qsuOld, 0.01);  // ~0.5° rotation
    item->setqsu(qsuPerturbed);
}
```

**Advantages**:
- Minimal code change
- Directly addresses the observed failure

**Disadvantages**:
- Perturbation magnitude is arbitrary
- May require multiple perturbations to escape singularity
- Doesn't address root cause

---

## 6. Recommendations

### 6.1 Short-Term Fix (Option B + C)

1. Before removing a constraint as redundant, check if it has a large violation
2. If violation > threshold, it's a configuration singularity
3. Apply perturbation and retry instead of removing the constraint

### 6.2 Long-Term Fix (Option A)

Implement Levenberg-Marquardt damping in the Newton-Raphson solver:
1. Add damping parameter λ to the normal equations
2. Implement adaptive λ control based on objective function reduction
3. Remove the need to detect/handle singularities specially

---

## 7. Assumptions Required for Proposed Fixes

### For Option B (Residual Check at Detection)

- **Assumption**: Constraint values are computed and accessible at singularity detection time
- **Verification**: Yes, `aG` is updated after each `calcPostDynCorrectorIteration()`

### For Option C (Perturbation)

- **Assumption**: Small perturbations from 180° will move to non-singular configurations
- **Verification**: Yes, the Jacobian becomes full-rank for any angle ≠ 180°
- **Assumption**: Perturbation won't push toward wrong solution
- **Caution**: Quaternion constraints have unique solution, so this should be safe

### For Option A (Levenberg-Marquardt)

- **Assumption**: The problem is fundamentally a least-squares optimization
- **Verification**: Yes, Newton-Raphson for constraints is equivalent to minimizing ||C(q)||²
- **Assumption**: Adaptive λ will converge
- **Verification**: Standard LM theory guarantees convergence for smooth objectives

---

## 8. References

1. **Levenberg-Marquardt Algorithm**: [Wikipedia](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm)
   - "A positive regularization or damping parameter is specifically added to remedy the possibility of the Jacobian matrix becoming rank-deficient."

2. **Configuration Singularities in Multibody Systems**: [ANSYS Motion Theory](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/motion_theory/motion_theory_fundamentals_constraints.html)
   - "Singular configuration occurs when the constraint equations lead to a Jacobian matrix that is linearly dependent."

3. **Quaternion-Based Constraints**: [RBDL Joint Modeling](https://rbdl.github.io/df/dbe/joint_description.html)
   - "When using quaternion joint type the model does not suffer from singularities."

4. **Jacobian Rank Analysis**: [Academia Paper on Quaternion Jacobians](https://www.academia.edu/2458960/Application_of_quaternion_algebra_to_the_efficient_computation_of_jacobians_for_holonomic_rheonomic_constraints)
   - Discusses efficient computation of constraint Jacobians using quaternion algebra.

5. **ACM Paper on Quaternion RBD**: [ACM Transactions on Graphics](https://dl.acm.org/doi/10.1145/3730872)
   - "A Versatile Quaternion-Based Constrained Rigid Body Dynamics" - guarantees satisfaction of kinematic constraints.

---

## 9. Conclusion

The 180° anti-parallel failure is a **configuration singularity**, not a constraint redundancy. The current solver incorrectly conflates these two conditions because it uses pivot magnitude as a proxy for redundancy.

The mathematically correct solution is Levenberg-Marquardt damping, which regularizes the Jacobian and prevents singularity-related failures. A shorter-term fix can check residuals at detection time to distinguish singularity from redundancy.

The user's intuition is correct: quaternion constraints should have exactly one solution, and a two-part assembly should always be solvable. The bug is in the solver's handling of singular configurations, not in the constraint formulation.

---

## 10. Implementation Status

**Status**: ✅ Implemented and working (December 2025)

### 10.1 Approach Taken

Rather than Options A-C described above, we implemented **Option D: Deterministic 180° Correction**. This approach exploits the specific geometric structure of the 180° singularity to apply an exact correction, rather than using generic numerical techniques.

**Key insight**: If the relative quaternion is approximately a 180° rotation about axis **n̂**, applying the same 180° rotation to one of the parts transforms the relative quaternion to identity. This is proven mathematically in the companion document [mathematical-proof-180-degree-correction.md](mathematical-proof-180-degree-correction.md).

### 10.2 Implementation Location

The 180° detection and correction logic is implemented in:
- `PosICNewtonRaphson.cpp` - Main detection/correction in `handleSingularMatrix()`
- `PosICNewtonRaphson.h` - Tracking data structures

### 10.3 What Works

1. **Detection**: Identifies 180° singularity when:
   - Relative quaternion is pure imaginary: `|qRel.w| < 0.1`
   - Constraint has large residual: `|residual| > 0.1`

2. **Correction**: Applies 180° rotation about dominant axis to part J

3. **Marker transforms**: Correction is conjugated by marker quaternion to transform from end-frame to part-frame: `qCorr_adj = qMarker × qCorr × conj(qMarker)`

4. **Shared parts**: Tracks corrected parts via `std::set<PartFrame*>` to prevent infinite loops when multiple joints share the same part

5. **Mini-solves**: Two-part assemblies with 180° anti-parallel orientation now solve correctly

### 10.4 Remaining Limitations

The 180° correction handles the configuration singularity, but cannot fix:
- **Geometrically impossible constraints**: If the assembly has conflicting constraints (e.g., over-constrained geometry), the solver will still report INCONSISTENT_CONSTRAINTS after exhausting retries
- **True redundancy**: Constraints that are algebraically redundant (zero residual) are still handled by the existing redundancy removal logic

### 10.5 References

See [mathematical-proof-180-degree-correction.md](mathematical-proof-180-degree-correction.md) for:
- Complete mathematical proof of the correction method
- Analysis of multiple singularity propagation
- Implementation algorithm details
