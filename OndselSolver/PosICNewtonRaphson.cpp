/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <assert.h>
#include <exception>
#include <map>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "PosICNewtonRaphson.h"
#include "SingularMatrixError.h"
#include "InconsistentConstraintsError.h"
#include "SystemSolver.h"
#include "Part.h"
#include "Constraint.h"
#include "CREATE.h"
#include "GESpMatParPvPrecise.h"
#include "GESpMatFullPvPosIC.h"
#include "Joint.h"

using namespace MbD;

void PosICNewtonRaphson::run()
{
	while (true) {
		try {
			//VectorNewtonRaphson::run();   //Inline to help debugging
			preRun();
			initializeLocally();
			initializeGlobally();
			iterate();
			postRun();
			break;
		}
		catch (InconsistentConstraintsError& ex) {
			auto inconsistentEqnNos = ex.getInconsistentEqnNos();
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->reactivateRedundantConstraints(); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(qsuOld); });

			// Build diagnostic information
			auto diagnostic = std::make_shared<InconsistencyDiagnostic>();
			std::map<std::string, int> partOccurrences;

			// Collect joint diagnostics
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
				auto joint = std::dynamic_pointer_cast<Joint>(item);
				if (joint) {
					auto jointDiag = joint->getJointDiagnostic(inconsistentEqnNos);
					if (!jointDiag.inconsistentConstraints.empty()) {
						diagnostic->joints.push_back(jointDiag);
						if (!jointDiag.partIName.empty()) partOccurrences[jointDiag.partIName]++;
						if (!jointDiag.partJName.empty()) partOccurrences[jointDiag.partJName]++;
					}
				}
			});

			// Identify affected part (appears most frequently)
			int maxCount = 0;
			for (const auto& pair : partOccurrences) {
				if (pair.second > maxCount) {
					maxCount = pair.second;
					diagnostic->affectedPartName = pair.first;
				}
			}

			// Compute total violation
			for (const auto& jointDiag : diagnostic->joints) {
				for (const auto& conDiag : jointDiag.inconsistentConstraints) {
					diagnostic->totalViolation += std::abs(conDiag.violation);
				}
			}

			// Build diagnostic YAML message with BEGIN/END markers for Python parsing
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(6);
			oss << "---BEGIN:INCONSISTENT_CONSTRAINTS---\n";
			oss << "Constraints are geometrically inconsistent (no solution exists)\n";
			oss << "  affected_part: \"" << diagnostic->affectedPartName << "\"\n";
			oss << "  total_violation: " << diagnostic->totalViolation << "\n";
			oss << "  joints:\n";
			for (const auto& jointDiag : diagnostic->joints) {
				oss << "    - name: \"" << jointDiag.name << "\"\n";
				oss << "      type: \"" << jointDiag.type << "\"\n";
				oss << "      part_i: \"" << jointDiag.partIName << "\"\n";
				oss << "      part_j: \"" << jointDiag.partJName << "\"\n";
				oss << "      lcs_i:\n";
				oss << "        name: \"" << jointDiag.lcsI.name << "\"\n";
				oss << "        position_on_part: [" << jointDiag.lcsI.positionOnPart[0] << ", "
					<< jointDiag.lcsI.positionOnPart[1] << ", " << jointDiag.lcsI.positionOnPart[2] << "]\n";
				oss << "        world_position: [" << jointDiag.lcsI.worldPosition[0] << ", "
					<< jointDiag.lcsI.worldPosition[1] << ", " << jointDiag.lcsI.worldPosition[2] << "]\n";
				oss << "      lcs_j:\n";
				oss << "        name: \"" << jointDiag.lcsJ.name << "\"\n";
				oss << "        position_on_part: [" << jointDiag.lcsJ.positionOnPart[0] << ", "
					<< jointDiag.lcsJ.positionOnPart[1] << ", " << jointDiag.lcsJ.positionOnPart[2] << "]\n";
				oss << "        world_position: [" << jointDiag.lcsJ.worldPosition[0] << ", "
					<< jointDiag.lcsJ.worldPosition[1] << ", " << jointDiag.lcsJ.worldPosition[2] << "]\n";
				oss << "      inconsistent_constraints:\n";
				for (const auto& conDiag : jointDiag.inconsistentConstraints) {
					oss << "        - type: \"" << conDiag.type << "\"\n";
					oss << "          violation: " << conDiag.violation << "\n";
				}
			}
			oss << "---END:INCONSISTENT_CONSTRAINTS---\n";

			// Re-throw with detailed message
			throw InconsistentConstraintsError(oss.str());
		}
		catch (const SingularMatrixError& ex) {
			auto redundantEqnNos = ex.getRedundantEqnNos();
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->removeRedundantConstraints(redundantEqnNos); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->constraintsReport(); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(qsuOld); });
		}
	}
}

void PosICNewtonRaphson::preRun()
{
	std::string str("MbD: Assembling system. ");
	system->logString(str);
	PosNewtonRaphson::preRun();
}

void PosICNewtonRaphson::assignEquationNumbers()
{
	auto parts = system->parts();
	//auto contactEndFrames = system->contactEndFrames();
	//auto uHolders = system->uHolders();
	auto essentialConstraints = system->essentialConstraints();
	auto displacementConstraints = system->displacementConstraints();
	auto perpendicularConstraints = system->perpendicularConstraints();
	size_t eqnNo = 0;
	for (auto& part : *parts) {
		part->iqX(eqnNo);
		eqnNo = eqnNo + 3;
		part->iqE(eqnNo);
		eqnNo = eqnNo + 4;
	}
	//for (auto& endFrm : *contactEndFrames) {
	//	endFrm->is(eqnNo);
	//	eqnNo = eqnNo + endFrm->sSize();
	//}
	//for (auto& uHolder : *uHolders) {
	//	uHolder->iu(eqnNo);
	//	eqnNo += 1;
	//}
	auto nEqns = eqnNo;	//C++ uses index 0.
	nqsu = nEqns;
	for (auto& con : *essentialConstraints) {
		con->iG = eqnNo;
		eqnNo += 1;
	}
	auto lastEssenConEqnNo = eqnNo - 1;
	for (auto& con : *displacementConstraints) {
		con->iG = eqnNo;
		eqnNo += 1;
	}
	auto lastDispConEqnNo = eqnNo - 1;
	for (auto& con : *perpendicularConstraints) {
		con->iG = eqnNo;
		eqnNo += 1;
	}
	auto lastEqnNo = eqnNo - 1;
	nEqns = eqnNo;	//C++ uses index 0.
	n = nEqns;
	auto rangelimits = { lastEssenConEqnNo + 1, lastDispConEqnNo + 1, lastEqnNo + 1 };
	pivotRowLimits = std::make_shared<std::vector<size_t>>(rangelimits);
}

bool PosICNewtonRaphson::isConverged()
{
	return this->isConvergedToNumericalLimit();
}

void PosICNewtonRaphson::handleSingularMatrix()
{
	nSingularMatrixError++;
	if (nSingularMatrixError == 1) {
		this->lookForRedundantConstraints();
		matrixSolver = this->matrixSolverClassNew();
	}
	else {
        auto& r = *matrixSolver;
		std::string str = typeid(r).name();
		if (str.find("GESpMatParPvMarkoFast") != std::string::npos) {
		    matrixSolver = CREATE<GESpMatParPvPrecise>::With();
		    this->solveEquations();
		}
		else {
            auto& msRef = *matrixSolver.get(); // extrapolated to suppress warning
            str = typeid(msRef).name();
            (void) msRef;                      // also for warning suppression
			if (str.find("GESpMatParPvPrecise") != std::string::npos) {
				this->lookForRedundantConstraints();
				matrixSolver = this->matrixSolverClassNew();
			} else {
				throw SimulationStoppingError("To be implemented.");
			}
		}
	}
}

void PosICNewtonRaphson::lookForRedundantConstraints()
{
	std::string str("MbD: Checking for redundant constraints.");
	system->logString(str);
	auto posICsolver = CREATE<GESpMatFullPvPosIC>::With();
	posICsolver->system = this;
	dx = posICsolver->solvewithsaveOriginal(pypx, y->negated(), false);
}
