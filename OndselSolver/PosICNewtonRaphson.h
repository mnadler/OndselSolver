/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include "AnyPosICNewtonRaphson.h"

namespace MbD {
    class PosICNewtonRaphson : public AnyPosICNewtonRaphson
    {
      //IC with over, fully or under constrained system
      //Perform redundant constraint removal for over constrained system
      //pivotRowLimits
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
        // Flag to prevent infinite retry loop when perturbing for anti-parallel detection
        bool hasRetriedWithPerturbation = false;
    };
}

