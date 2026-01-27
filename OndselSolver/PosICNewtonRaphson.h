/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include <set>
#include <limits>
#include "AnyPosICNewtonRaphson.h"

namespace MbD {
    class Constraint;  // Forward declaration
    class PartFrame;   // Forward declaration

    class PosICNewtonRaphson : public AnyPosICNewtonRaphson
    {
      //IC with over, fully or under constrained system
      //Perform redundant constraint removal for over constrained system
      //Uses iterative constraint protection to distinguish truly redundant from essential constraints
      //pivotRowLimits protectedConstraints
    public:
        PosICNewtonRaphson(){}

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
        // Best effort converged state - saved before retrying, restored before throwing
        FColDsptr bestEffortState;
        // Track the best yNorm seen during iteration (lowest = best solution)
        double bestYNorm = std::numeric_limits<double>::max();
        // Track parts that have been corrected for 180° singularity (prevents flip-flopping)
        std::set<PartFrame*> correctedPartsFor180;

        // Clear correction tracking - call at start of each NEW solve, not between retries
        // See docs/fix-180-degree-correction-oscillation.md for rationale
        void clearCorrectedParts() { correctedPartsFor180.clear(); }
    };
}

