# Quaternion (Euler Parameter) Implementation in OndselSolver

## Overview

OndselSolver uses **Euler parameters** (unit quaternions) to represent 3D rotations. This approach avoids the singularities inherent in Euler angles (gimbal lock) and provides smooth, continuous rotation representation throughout the entire rotation space.

A quaternion **q = [e₀, e₁, e₂, e₃]** encodes a rotation where:
- **e₀, e₁, e₂** are the imaginary components (rotation axis × sin(θ/2))
- **e₃** is the scalar component (cos(θ/2))

## Why Quaternions Instead of Euler Angles?

| Aspect | Euler Angles | Quaternions |
|--------|--------------|-------------|
| Parameters | 3 angles (φ, θ, ψ) | 4 components (e₀, e₁, e₂, e₃) |
| Singularities | Gimbal lock at θ = ±90° | None |
| Interpolation | Discontinuous | Smooth (SLERP) |
| Composition | Matrix multiplication | Quaternion multiplication |
| Constraint | None | e₀² + e₁² + e₂² + e₃² = 1 |

The trade-off is one additional variable (4 vs 3), but the elimination of singularities is critical for robust numerical solving.

## Mathematical Formulation

### Quaternion Representation

A rotation of angle **θ** about unit axis **n = (nₓ, nᵧ, n_z)** is encoded as:

```
q = [e₀, e₁, e₂, e₃] = [nₓ·sin(θ/2), nᵧ·sin(θ/2), n_z·sin(θ/2), cos(θ/2)]
```

**Note:** OndselSolver uses 0-indexed arrays in C++, so `qE->at(0)` = e₀, etc.

### Unit Quaternion Constraint

Quaternions must remain normalized:

```
G = e₀² + e₁² + e₂² + e₃² - 1 = 0
```

This is enforced by `EulerConstraint` with Jacobian:

```
∂G/∂qE = [2e₀, 2e₁, 2e₂, 2e₃]
```

### Rotation Matrix from Quaternion

The 3×3 rotation matrix **A** (Direction Cosine Matrix) is computed as:

```
A = B · Cᵀ
```

Where **B** and **C** are 3×4 matrices:

```
B = | e₃  -e₂   e₁  -e₀ |
    | e₂   e₃  -e₀  -e₁ |
    |-e₁   e₀   e₃  -e₂ |

C = | e₃   e₂  -e₁  -e₀ |
    |-e₂   e₃   e₀  -e₁ |
    | e₁  -e₀   e₃  -e₂ |
```

The resulting rotation matrix:

```
A = | 1-2(e₁²+e₂²)    2(e₀e₁-e₂e₃)    2(e₀e₂+e₁e₃) |
    | 2(e₀e₁+e₂e₃)    1-2(e₀²+e₂²)    2(e₁e₂-e₀e₃) |
    | 2(e₀e₂-e₁e₃)    2(e₁e₂+e₀e₃)    1-2(e₀²+e₁²) |
```

### Quaternion to Rotation Matrix (Implementation)

```cpp
// In EulerParameters<T>::calcABC()
void calcABC() {
    // Build aB matrix (3x4)
    aB = std::make_shared<FullMatrix<double>>(3, 4);
    aB->at(0, 0) =  e3; aB->at(0, 1) = -e2; aB->at(0, 2) =  e1; aB->at(0, 3) = -e0;
    aB->at(1, 0) =  e2; aB->at(1, 1) =  e3; aB->at(1, 2) = -e0; aB->at(1, 3) = -e1;
    aB->at(2, 0) = -e1; aB->at(2, 1) =  e0; aB->at(2, 2) =  e3; aB->at(2, 3) = -e2;

    // Build aC matrix (3x4)
    aC = std::make_shared<FullMatrix<double>>(3, 4);
    aC->at(0, 0) =  e3; aC->at(0, 1) =  e2; aC->at(0, 2) = -e1; aC->at(0, 3) = -e0;
    aC->at(1, 0) = -e2; aC->at(1, 1) =  e3; aC->at(1, 2) =  e0; aC->at(1, 3) = -e1;
    aC->at(2, 0) =  e1; aC->at(2, 1) = -e0; aC->at(2, 2) =  e3; aC->at(2, 3) = -e2;

    // Compute rotation matrix: A = B · Cᵀ
    aA = aB->timesTransposeFullMatrix(aC);
}
```

### Rotation Matrix to Quaternion (Shepperd's Method)

Converting a rotation matrix back to a quaternion requires care to avoid numerical instability. OndselSolver uses **Shepperd's method**, which selects the computation path based on which quaternion component is largest:

```cpp
// EulerParameters<double>::fromRotationMatrix() - static factory method
static std::shared_ptr<EulerParameters<double>> fromRotationMatrix(FMatDsptr aA) {
    auto qE = std::make_shared<EulerParameters<double>>(4);
    double trace = aA->at(0)->at(0) + aA->at(1)->at(1) + aA->at(2)->at(2);

    if (trace > 0) {
        // e₃ (scalar) is largest - compute it first
        double s = 0.5 / std::sqrt(trace + 1.0);
        qE->at(3) = 0.25 / s;  // w
        qE->at(0) = (aA->at(2)->at(1) - aA->at(1)->at(2)) * s;  // x
        qE->at(1) = (aA->at(0)->at(2) - aA->at(2)->at(0)) * s;  // y
        qE->at(2) = (aA->at(1)->at(0) - aA->at(0)->at(1)) * s;  // z
    }
    else if (aA->at(0)->at(0) > aA->at(1)->at(1) && aA->at(0)->at(0) > aA->at(2)->at(2)) {
        // e₀ (x) is largest
        double s = 2.0 * std::sqrt(1.0 + aA->at(0)->at(0) - aA->at(1)->at(1) - aA->at(2)->at(2));
        qE->at(3) = (aA->at(2)->at(1) - aA->at(1)->at(2)) / s;
        qE->at(0) = 0.25 * s;
        qE->at(1) = (aA->at(0)->at(1) + aA->at(1)->at(0)) / s;
        qE->at(2) = (aA->at(0)->at(2) + aA->at(2)->at(0)) / s;
    }
    else if (aA->at(1)->at(1) > aA->at(2)->at(2)) {
        // e₁ (y) is largest
        double s = 2.0 * std::sqrt(1.0 + aA->at(1)->at(1) - aA->at(0)->at(0) - aA->at(2)->at(2));
        qE->at(3) = (aA->at(0)->at(2) - aA->at(2)->at(0)) / s;
        qE->at(0) = (aA->at(0)->at(1) + aA->at(1)->at(0)) / s;
        qE->at(1) = 0.25 * s;
        qE->at(2) = (aA->at(1)->at(2) + aA->at(2)->at(1)) / s;
    }
    else {
        // e₂ (z) is largest
        double s = 2.0 * std::sqrt(1.0 + aA->at(2)->at(2) - aA->at(0)->at(0) - aA->at(1)->at(1));
        qE->at(3) = (aA->at(1)->at(0) - aA->at(0)->at(1)) / s;
        qE->at(0) = (aA->at(0)->at(2) + aA->at(2)->at(0)) / s;
        qE->at(1) = (aA->at(1)->at(2) + aA->at(2)->at(1)) / s;
        qE->at(2) = 0.25 * s;
    }
    return qE;
}
```

This avoids division by small numbers and maintains numerical precision.

## Angular Velocity and Quaternion Rates

### From Angular Velocity to Quaternion Rate

Given angular velocity **ω** in body frame:

```
q̇ = ½ · Bᵀ · ω
```

Implementation:
```cpp
// EulerParametersDot::FromqEOpAndOmegaOpO()
qEdot = aB->transposeTimesFullColumn(omega->times(0.5));
```

### From Quaternion Rate to Angular Velocity

```
ω = 2 · B · q̇
```

Implementation:
```cpp
// EulerParametersDot::omeOpO()
FColDsptr omega = aB()->timesFullColumn(this)->times(2.0);
```

### Quaternion Rate Constraint

The quaternion rate must be perpendicular to the quaternion itself:

```
q · q̇ = 0
```

This is automatically satisfied when computing q̇ from ω using the formula above.

## World Orientation Quaternion (qEO)

### The Problem: Part vs. End Frame Orientation

A joint connects two **end frames** (marker frames), not parts directly. Each end frame has:
- A **part quaternion** `qE_part` - the part's orientation in world coordinates
- A **marker quaternion** `qE_marker` - the marker's orientation relative to the part

The **world orientation** of the end frame is the composition:

```
qE_world = qE_part * qE_marker
```

### Implementation: EndFrameqc::qEO()

```cpp
std::shared_ptr<EulerParameters<double>> EndFrameqc::qEO()
{
    // Compute world orientation quaternion: qE_world = qE_part * qE_marker
    auto qEpart = markerFrame->qE();
    auto qEmarker = markerFrame->qEpm;

    auto result = std::make_shared<EulerParameters<double>>(4);

    // Hamilton product: p = q1 * q2
    // OndselSolver convention: [e0, e1, e2, e3] = [x, y, z, w]
    double q1_0 = qEpart->at(0), q1_1 = qEpart->at(1), q1_2 = qEpart->at(2), q1_3 = qEpart->at(3);
    double q2_0 = qEmarker->at(0), q2_1 = qEmarker->at(1), q2_2 = qEmarker->at(2), q2_3 = qEmarker->at(3);

    result->at(0) = q1_3*q2_0 + q1_0*q2_3 + q1_1*q2_2 - q1_2*q2_1;  // x
    result->at(1) = q1_3*q2_1 - q1_0*q2_2 + q1_1*q2_3 + q1_2*q2_0;  // y
    result->at(2) = q1_3*q2_2 + q1_0*q2_1 - q1_1*q2_0 + q1_2*q2_3;  // z
    result->at(3) = q1_3*q2_3 - q1_0*q2_0 - q1_1*q2_1 - q1_2*q2_2;  // w

    return result;
}
```

This is critical for quaternion constraints to work correctly - they must compare the actual end frame orientations, not just the part orientations.

## Orientation Constraints: Direction Cosine vs. Quaternion

### The Anti-Parallel Problem with Direction Cosine Constraints

Traditional MBD solvers use **direction cosine constraints** to enforce relative orientation between frames. For a fixed joint (all axes aligned), three perpendicularity constraints are used:

```cpp
// OLD approach: Direction cosine perpendicularity constraints
// Constraint: axis_i of frame I is perpendicular to axis_j of frame J
DirectionCosineConstraintIqcJqc(frmI, frmJ, 1, 0);  // Y_I ⊥ X_J
DirectionCosineConstraintIqcJqc(frmI, frmJ, 2, 0);  // Z_I ⊥ X_J
DirectionCosineConstraintIqcJqc(frmI, frmJ, 2, 1);  // Z_I ⊥ Y_J
```

**The problem:** Perpendicularity constraints have **TWO valid solutions**:
1. Frames aligned (0° rotation) ✓
2. Frames anti-parallel (180° rotation) ✗

When parts start in an "exploded" configuration (all at origin with identity orientation), the solver may converge to the **wrong** solution where parts are 180° rotated from their correct orientation.

### The Solution: Quaternion-Based Constraints

Quaternion constraints directly constrain the **relative quaternion** between frames:

```
q_rel = conj(q_I) * q_J
```

For aligned frames, q_rel should be [0, 0, 0, 1] (identity rotation). The constraint enforces that the **imaginary components are zero**:

```
Im(q_rel) = [0, 0, 0]
```

This has only **ONE valid solution**: aligned orientation.

### Why Quaternion Constraints Work

| Constraint Type | Valid Solutions | Anti-parallel? |
|-----------------|-----------------|----------------|
| Direction cosine perpendicularity | 2 (aligned, anti-parallel) | Accepted ✗ |
| Quaternion imaginary = 0 | 1 (aligned only) | Rejected ✓ |

The mathematical reason:
- **Direction cosines** are elements of the rotation matrix A
- A is the **same** for rotations q and -q (180° apart)
- **Quaternion constraints** distinguish between q and -q
- conj(q_I) * q_J = [0,0,0,1] only when frames are truly aligned

### Implementation in FixedJoint

```cpp
void MbD::FixedJoint::initializeGlobally()
{
    if (constraints->empty())
    {
        createAtPointConstraints();  // Position constraints (unchanged)

        // NEW: Quaternion-based orientation constraints
        // Each constrains one imaginary component of conj(qI) * qJ to zero
        addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 0));  // x = 0
        addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 1));  // y = 0
        addConstraint(CREATE<QuaternionConstraintIqcJqc>::With(frmI, frmJ, 2));  // z = 0

        this->root()->hasChanged = true;
    }
    // ...
}
```

## QuaternionConstraint Class Hierarchy

### Overview

The quaternion constraint classes mirror the direction cosine constraint hierarchy:

```
ConstraintIJ
  └── QuaternionConstraintIJ          // Base: constrains Im[axis](conj(qI)*qJ) = 0
       ├── QuaternionConstraintIqcJc   // Frame I has qc (rotating), J is constant
       └── QuaternionConstraintIqcJqc  // Both frames have qc (both rotating)

KinematicIeJe
  └── QuaternionIecJec                 // Computes Im[axis](conj(qI)*qJ) value
       ├── QuaternionIeqcJec           // With Jacobians for rotating I
       └── QuaternionIeqcJeqc          // With Jacobians for both rotating
```

### Naming Convention

- **I/J**: Frame I and Frame J
- **e**: End frame
- **c**: Constant (not differentiated)
- **qc**: Has quaternion coordinates (rotating, needs Jacobians)

Examples:
- `QuaternionIecJec`: Both frames constant (e.g., ground-to-ground)
- `QuaternionIeqcJec`: Frame I rotating, Frame J constant
- `QuaternionIeqcJeqc`: Both frames rotating (most common)

### QuaternionIecJec: The Core Computation

```cpp
void QuaternionIecJec::calcPostDynCorrectorIteration()
{
    // Get world orientation quaternions (accounts for marker frame)
    auto efrmI = std::static_pointer_cast<EndFrameqc>(frmI);
    auto efrmJ = std::static_pointer_cast<EndFrameqc>(frmJ);
    auto qEI = efrmI->qEO();  // World quaternion of frame I
    auto qEJ = efrmJ->qEO();  // World quaternion of frame J

    double qI0 = qEI->at(0), qI1 = qEI->at(1), qI2 = qEI->at(2), qI3 = qEI->at(3);
    double qJ0 = qEJ->at(0), qJ1 = qEJ->at(1), qJ2 = qEJ->at(2), qJ3 = qEJ->at(3);

    // Hamilton product: p = conj(qI) * qJ
    // conj(qI) = (-qI0, -qI1, -qI2, qI3)
    // p[0] = qI3*qJ0 - qI0*qJ3 - qI1*qJ2 + qI2*qJ1  (x imaginary)
    // p[1] = qI3*qJ1 + qI0*qJ2 - qI1*qJ3 - qI2*qJ0  (y imaginary)
    // p[2] = qI3*qJ2 - qI0*qJ1 + qI1*qJ0 - qI2*qJ3  (z imaginary)

    if (axis == 0) {
        aQijIeJe = qI3*qJ0 - qI0*qJ3 - qI1*qJ2 + qI2*qJ1;
    }
    else if (axis == 1) {
        aQijIeJe = qI3*qJ1 + qI0*qJ2 - qI1*qJ3 - qI2*qJ0;
    }
    else { // axis == 2
        aQijIeJe = qI3*qJ2 - qI0*qJ1 + qI1*qJ0 - qI2*qJ3;
    }
}
```

### QuaternionConstraintIqcJqc: With Jacobians

For rotating frames, the solver needs Jacobians ∂G/∂qE. The `QuaternionIeqcJeqc` class computes:
- `pAQIeJepEI`: ∂(Im[axis])/∂qE_I (4 values)
- `pAQIeJepEJ`: ∂(Im[axis])/∂qE_J (4 values)

These are assembled into the system Jacobian for Newton-Raphson iteration.

## Class Hierarchy (Complete)

```
FullVector<T>
  └── EulerArray<T>               // Base for quaternion arrays
       ├── EulerParameters<T>         // Quaternion q = [e₀, e₁, e₂, e₃]
       │    - aA: rotation matrix (3×3)
       │    - aB, aC: helper matrices (3×4)
       │    - pApE: ∂A/∂qE derivatives
       │    + fromRotationMatrix(): Shepperd's method
       │
       ├── EulerParametersDot<T>      // Quaternion velocity q̇
       │    - aBdot, aCdot: time derivatives of B, C
       │    - omeOpO(): angular velocity extraction
       │
       └── EulerParametersDDot<T>     // Quaternion acceleration q̈
            - alpOpO(): angular acceleration extraction

Constraint
  ├── EulerConstraint                 // q·q = 1 normalization
  ├── DirectionCosineConstraintIJ     // Rotation matrix element constraints (OLD)
  └── QuaternionConstraintIJ          // Relative quaternion constraints (NEW)
       ├── QuaternionConstraintIqcJc
       └── QuaternionConstraintIqcJqc

KinematicIeJe
  ├── DirectionCosineIecJec           // Computes A_I column · A_J column
  └── QuaternionIecJec                // Computes Im[axis](conj(qI)*qJ)
       ├── QuaternionIeqcJec
       └── QuaternionIeqcJeqc
```

## PartFrame: Quaternion State Storage

Each rigid body has a `PartFrame` that stores its kinematic state:

```cpp
class PartFrame {
    // Position (translation)
    FColDsptr qX;           // [x, y, z]
    FColDsptr qXdot;        // [ẋ, ẏ, ż]
    FColDsptr qXddot;       // [ẍ, ÿ, z̈]

    // Orientation (quaternion)
    std::shared_ptr<EulerParameters<double>> qE;      // [e₀, e₁, e₂, e₃]
    std::shared_ptr<EulerParametersDot<double>> qEdot; // [ė₀, ė₁, ė₂, ė₃]
    FColDsptr qEddot;                                  // [ë₀, ë₁, ë₂, ë₃]

    // Constraints
    std::shared_ptr<Constraint> aGeu;  // Euler parameter constraint (normalization)
    std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> aGabs;  // Additional constraints

    // System variable indices
    size_t iqX;  // Index of qX in system state vector
    size_t iqE;  // Index of qE in system state vector
};
```

## Handling Anti-Parallel Configurations (180° Rotations)

### Two-Pronged Approach

OndselSolver addresses the anti-parallel problem at two levels:

**1. Quaternion Constraints (Primary Fix)**
- Replace direction cosine constraints with quaternion constraints
- Eliminates anti-parallel solutions mathematically
- Used in `FixedJoint` (and applicable to other joints)

**2. Iterative Constraint Protection (Backup)**
- Handles cases where Jacobian becomes singular at 180°
- Learns which constraints are essential vs. truly redundant
- Protects essential constraints from removal during retry

See [plan-iterative-constraint-protection.md](plan-iterative-constraint-protection.md) for full details on the iterative approach.

### When Each Approach Applies

| Scenario | Primary Fix | Backup |
|----------|-------------|--------|
| FixedJoint with anti-parallel start | Quaternion constraints | — |
| RevoluteJoint at 180° rotation | Direction cosine (still used) | Iterative protection |
| Inconsistent assembly | — | Iterative protection + best-effort |

## Derivatives and Jacobians

### Partial Derivative of Rotation Matrix w.r.t. Quaternion

The solver needs **∂A/∂qE** for computing constraint Jacobians. This is a 3×3×4 tensor, but is typically accessed column-by-column:

```cpp
// In EulerParameters<T>::calcpApE()
// Computes ∂Aⱼ/∂qEᵢ for each column j and component i
```

### Chain Rule for Constraint Jacobians

For a constraint G(A(qE)):

```
∂G/∂qE = ∂G/∂A · ∂A/∂qE
```

For quaternion constraints G(qE_I, qE_J):

```
∂G/∂qE_I = direct differentiation of Hamilton product
∂G/∂qE_J = direct differentiation of Hamilton product
```

This is computed in `QuaternionIeqcJeqc`.

## Numerical Considerations

### Quaternion Conditioning

```cpp
// In FullVector<T>::conditionSelf()
void conditionSelf() {
    for (size_t i = 0; i < size(); i++) {
        if (std::abs(at(i)) < std::numeric_limits<double>::epsilon()) {
            at(i) = 0.0;
        }
    }
}
```

### Quaternion Normalization

```cpp
// In FullVector<T>::normalizeSelf()
void normalizeSelf() {
    double mag = length();
    if (mag > 0) {
        for (size_t i = 0; i < size(); i++) {
            at(i) /= mag;
        }
    }
}
```

### Avoiding Antipodal Ambiguity

Quaternions q and -q represent the same rotation. To maintain consistency:

```cpp
// Ensure scalar component is positive (canonical form)
if (qE->at(3) < 0) {
    qE->negateSelf();
}
```

## Summary

| Concept | Formula/Implementation |
|---------|----------------------|
| Quaternion | q = [e₀, e₁, e₂, e₃] with e₀²+e₁²+e₂²+e₃² = 1 |
| Rotation matrix | A = B · Cᵀ |
| Angular velocity from q̇ | ω = 2 · B · q̇ |
| q̇ from angular velocity | q̇ = ½ · Bᵀ · ω |
| World quaternion | qE_world = qE_part × qE_marker |
| Relative quaternion | q_rel = conj(q_I) × q_J |
| Euler constraint | ∂G/∂qE = [2e₀, 2e₁, 2e₂, 2e₃] |
| Quaternion constraint | Im[axis](conj(q_I) × q_J) = 0 |
| Matrix-to-quaternion | Shepperd's method (numerically stable) |
| Anti-parallel fix | Quaternion constraints (single solution) |

## Key Files

| File | Purpose |
|------|---------|
| `EulerParameters.h/cpp` | Quaternion class with rotation matrix computation |
| `EulerParametersDot.h/cpp` | Quaternion velocity, angular velocity conversion |
| `EulerParametersDDot.h/cpp` | Quaternion acceleration |
| `EulerConstraint.h/cpp` | Normalization constraint q·q = 1 |
| `PartFrame.h/cpp` | Stores quaternion state per rigid body |
| `EndFrameqc.h/cpp` | End frame with qEO() world quaternion method |
| `DirectionCosineConstraintIJ.h/cpp` | Joint constraints via rotation matrix elements (old) |
| `QuaternionConstraintIJ.h/cpp` | Joint constraints via relative quaternion (new) |
| `QuaternionIecJec.h/cpp` | Computes relative quaternion imaginary component |
| `QuaternionIeqcJeqc.h/cpp` | With Jacobians for both rotating frames |
| `FixedJoint.cpp` | Uses QuaternionConstraint instead of DirectionCosine |
| `PosICNewtonRaphson.cpp` | Iterative constraint protection for singularities |

## References

- Shepperd, S.W. (1978). "Quaternion from Rotation Matrix". *Journal of Guidance and Control*.
- Nikravesh, P.E. (1988). *Computer-Aided Analysis of Mechanical Systems*. Prentice-Hall.
- Haug, E.J. (1989). *Computer Aided Kinematics and Dynamics of Mechanical Systems*. Allyn and Bacon.
