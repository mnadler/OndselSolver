/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>

#include "GESpMatFullPvPosIC.h"
#include "SingularMatrixError.h"
#include "InconsistentConstraintsError.h"
#include "PosICNewtonRaphson.h"

using namespace MbD;

void GESpMatFullPvPosIC::preSolvewithsaveOriginal(SpMatDsptr spMat, FColDsptr fullCol, bool saveOriginal)
{
	GESpMatFullPv::preSolvewithsaveOriginal(spMat, fullCol, saveOriginal);
	if (system == nullptr) {
		pivotRowLimits = std::make_shared<std::vector<size_t>>();
	}
	else {
		pivotRowLimits = system->pivotRowLimits;
	}
	pivotRowLimit = SIZE_MAX;
}

void GESpMatFullPvPosIC::doPivoting(size_t p)
{
	//"Used by Gauss Elimination only."
	//"Swap rows but keep columns in place."
	//"The elements below the diagonal are removed column by column."

	double max = 0.0;
	auto pivotRow = p;
	auto pivotCol = p;
	for (size_t j = p; j < n; j++)
	{
		rowPositionsOfNonZerosInColumns->at(colOrder->at(j))->clear();
	}
	if (pivotRowLimit == SIZE_MAX || p >= pivotRowLimit) {
		pivotRowLimit = *std::find_if(
			pivotRowLimits->begin(), pivotRowLimits->end(),
			[&](size_t limit) { return limit > p; });
	}
	for (size_t i = p; i < pivotRowLimit; i++)
	{
		auto& rowi = matrixA->at(i);
		for (auto const& kv : *rowi) {
			rowPositionsOfNonZerosInColumns->at(kv.first)->push_back(i);
			auto aij = kv.second;
			auto mag = aij;
			if (mag < 0.0) mag = -mag;
			if (max < mag) {
				max = mag;
				pivotRow = i;
				pivotCol = positionsOfOriginalCols->at(kv.first);
			}
		}
	}
	if (p != pivotRow) {
		matrixA->swapElems(p, pivotRow);
		rightHandSideB->swapElems(p, pivotRow);
		rowOrder->swapElems(p, pivotRow);
	}
	if (p != pivotCol) {
		colOrder->swapElems(p, pivotCol);
		positionsOfOriginalCols->at(colOrder->at(p)) = p;
		positionsOfOriginalCols->at(colOrder->at(pivotCol)) = pivotCol;
	}
	pivotValues->at(p) = max;
	if (max < singularPivotTolerance) {
		auto itr = std::find_if(
			pivotRowLimits->begin(), pivotRowLimits->end(),
			[&](size_t limit) { return limit > pivotRowLimit; });
		if (itr == pivotRowLimits->end()) {
			auto begin = rowOrder->begin() + p;
			auto end = rowOrder->begin() + pivotRowLimit;
			auto eqnNos = std::make_shared<FullColumn<size_t>>(begin, end);

			// Check if the RHS entries for dependent rows are near-zero.
			// If not, the constraints are inconsistent (no solution exists),
			// not merely redundant (infinite solutions exist).
			// This implements the augmented matrix rank test: rank(A) < rank([A|b])
			// means inconsistent, while rank(A) = rank([A|b]) means redundant.
			bool isInconsistent = false;
			double consistencyTolerance = 1.0e-6;
			for (size_t i = p; i < pivotRowLimit; i++) {
				if (std::abs(rightHandSideB->at(rowOrder->at(i))) > consistencyTolerance) {
					isInconsistent = true;
					break;
				}
			}

			if (isInconsistent) {
				// Only include equations with non-zero RHS (truly inconsistent)
				auto inconsistentEqnNos = std::make_shared<FullColumn<size_t>>();
				auto rhsValues = std::make_shared<std::vector<double>>();
				for (size_t i = p; i < pivotRowLimit; i++) {
					double rhs = rightHandSideB->at(rowOrder->at(i));
					if (std::abs(rhs) > consistencyTolerance) {
						inconsistentEqnNos->push_back(rowOrder->at(i));
						rhsValues->push_back(rhs);
					}
				}
				throw InconsistentConstraintsError(
					"Constraints are geometrically inconsistent (no solution exists)", inconsistentEqnNos, rhsValues);
			}

			throwSingularMatrixError("", eqnNos);
		}
		else {
			pivotRowLimit = *itr;
		}
		return this->doPivoting(p);
	}
	auto jp = colOrder->at(p);
	rowPositionsOfNonZerosInPivotColumn = rowPositionsOfNonZerosInColumns->at(jp);
	for (size_t i = pivotRowLimit; i < m; i++)
	{
		auto& spRowi = matrixA->at(i);
		if (spRowi->find(jp) != spRowi->end()) {
			rowPositionsOfNonZerosInPivotColumn->push_back(i);
		}
	}
	if (rowPositionsOfNonZerosInPivotColumn->front() == p) {
		rowPositionsOfNonZerosInPivotColumn->erase(rowPositionsOfNonZerosInPivotColumn->begin());
	}
	markowitzPivotColCount = rowPositionsOfNonZerosInPivotColumn->size();
}