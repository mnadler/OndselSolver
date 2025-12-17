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
#include <set>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "PosICNewtonRaphson.h"
#include "SingularMatrixError.h"
#include "InconsistentConstraintsError.h"
#include "SimulationStoppingError.h"
#include "SystemSolver.h"
#include "Part.h"
#include "PartFrame.h"
#include "Constraint.h"
#include "RedundantConstraint.h"
#include "CREATE.h"
#include "GESpMatParPvPrecise.h"
#include "GESpMatFullPvPosIC.h"
#include "Joint.h"
#include "ConstraintIJ.h"
#include "EndFramec.h"

using namespace MbD;

void PosICNewtonRaphson::run()
{
	// Clear any previous tracking
	removedEqnNos = nullptr;
	removedRhsAtDetection = nullptr;

	while (true) {
		try {
			//VectorNewtonRaphson::run();   //Inline to help debugging
			preRun();
			initializeLocally();
			initializeGlobally();
			iterate();
			postRun();
			// After successful convergence, verify removed constraints are satisfied
			verifyRemovedConstraintsAtConvergence();
			break;
		}
		catch (InconsistentConstraintsError& ex) {
			auto inconsistentEqnNos = ex.getInconsistentEqnNos();
			auto rhsValues = ex.getRhsValues();
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->reactivateRedundantConstraints(); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(qsuOld); });

			// Build a map from equation number to RHS value for quick lookup
			std::map<size_t, double> eqnToRhs;
			if (inconsistentEqnNos && rhsValues) {
				for (size_t i = 0; i < inconsistentEqnNos->size(); i++) {
					eqnToRhs[inconsistentEqnNos->at(i)] = rhsValues->at(i);
				}
			}

			// Build diagnostic information
			auto diagnostic = std::make_shared<InconsistencyDiagnostic>();
			std::map<std::string, int> partOccurrences;
			std::set<size_t> matchedEqnNos;  // Track which equation numbers were matched

			// Collect joint diagnostics
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
				auto joint = std::dynamic_pointer_cast<Joint>(item);
				if (joint) {
					auto jointDiag = joint->getJointDiagnostic(inconsistentEqnNos, rhsValues);
					if (!jointDiag.inconsistentConstraints.empty()) {
						diagnostic->joints.push_back(jointDiag);
						if (!jointDiag.partIName.empty()) partOccurrences[jointDiag.partIName]++;
						if (!jointDiag.partJName.empty()) partOccurrences[jointDiag.partJName]++;
						for (const auto& conDiag : jointDiag.inconsistentConstraints) {
							matchedEqnNos.insert(conDiag.equationNumber);
						}
					}
				}
			});

			// Collect Part constraint diagnostics (aGeu, aGabs)
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
				auto part = std::dynamic_pointer_cast<Part>(item);
				if (part && part->partFrame) {
					// Check Euler constraint (aGeu)
					auto aGeu = part->partFrame->aGeu;
					if (aGeu && inconsistentEqnNos) {
						auto it = std::find(inconsistentEqnNos->begin(), inconsistentEqnNos->end(), aGeu->iG);
						if (it != inconsistentEqnNos->end()) {
							PartConstraintDiagnostic pcDiag;
							pcDiag.partName = part->name;
							pcDiag.constraintType = aGeu->constraintSpec();
							pcDiag.equationNumber = aGeu->iG;
							// Use RHS value from exception if available, otherwise use current aG
							auto rhsIt = eqnToRhs.find(aGeu->iG);
							pcDiag.violation = (rhsIt != eqnToRhs.end()) ? rhsIt->second : aGeu->aG;
							diagnostic->partConstraints.push_back(pcDiag);
							partOccurrences[part->name]++;
							matchedEqnNos.insert(aGeu->iG);
						}
					}
					// Check absolute constraints (aGabs)
					if (part->partFrame->aGabs && inconsistentEqnNos) {
						for (auto& aGab : *(part->partFrame->aGabs)) {
							auto it = std::find(inconsistentEqnNos->begin(), inconsistentEqnNos->end(), aGab->iG);
							if (it != inconsistentEqnNos->end()) {
								PartConstraintDiagnostic pcDiag;
								pcDiag.partName = part->name;
								pcDiag.constraintType = aGab->constraintSpec();
								pcDiag.equationNumber = aGab->iG;
								auto rhsIt = eqnToRhs.find(aGab->iG);
								pcDiag.violation = (rhsIt != eqnToRhs.end()) ? rhsIt->second : aGab->aG;
								diagnostic->partConstraints.push_back(pcDiag);
								partOccurrences[part->name]++;
								matchedEqnNos.insert(aGab->iG);
							}
						}
					}
				}
			});

			// Collect unmatched equation numbers (ones that weren't found in any constraint)
			std::vector<std::pair<size_t, double>> unmatchedEqns;
			if (inconsistentEqnNos) {
				for (size_t i = 0; i < inconsistentEqnNos->size(); i++) {
					size_t eqnNo = inconsistentEqnNos->at(i);
					if (matchedEqnNos.find(eqnNo) == matchedEqnNos.end()) {
						double rhs = (rhsValues && i < rhsValues->size()) ? rhsValues->at(i) : 0.0;
						unmatchedEqns.push_back({eqnNo, rhs});
					}
				}
			}

			// Identify affected part (appears most frequently)
			int maxCount = 0;
			for (const auto& pair : partOccurrences) {
				if (pair.second > maxCount) {
					maxCount = pair.second;
					diagnostic->affectedPartName = pair.first;
				}
			}

			// Compute total violation from joints
			for (const auto& jointDiag : diagnostic->joints) {
				for (const auto& conDiag : jointDiag.inconsistentConstraints) {
					diagnostic->totalViolation += std::abs(conDiag.violation);
				}
			}
			// Compute total violation from part constraints
			for (const auto& pcDiag : diagnostic->partConstraints) {
				diagnostic->totalViolation += std::abs(pcDiag.violation);
			}
			// Compute total violation from unmatched equations
			for (const auto& unmatched : unmatchedEqns) {
				diagnostic->totalViolation += std::abs(unmatched.second);
			}

			// Build diagnostic YAML message with BEGIN/END markers for Python parsing
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(6);
			oss << "---BEGIN:INCONSISTENT_CONSTRAINTS---\n";
			oss << "Constraints are geometrically inconsistent (no solution exists)\n";
			oss << "  affected_part: \"" << diagnostic->affectedPartName << "\"\n";
			oss << "  total_violation: " << diagnostic->totalViolation << "\n";
			oss << "  inconsistent_equation_count: " << (inconsistentEqnNos ? inconsistentEqnNos->size() : 0) << "\n";
			oss << "  nqsu: " << nqsu << "\n";
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
			oss << "  part_constraints:\n";
			for (const auto& pcDiag : diagnostic->partConstraints) {
				oss << "    - part: \"" << pcDiag.partName << "\"\n";
				oss << "      type: \"" << pcDiag.constraintType << "\"\n";
				oss << "      equation_number: " << pcDiag.equationNumber << "\n";
				oss << "      violation: " << pcDiag.violation << "\n";
			}
			oss << "  unmatched_equations:\n";
			for (const auto& unmatched : unmatchedEqns) {
				oss << "    - equation_number: " << unmatched.first << "\n";
				oss << "      rhs_value: " << unmatched.second << "\n";
				// Provide hint about what this equation might be
				if (unmatched.first < nqsu) {
					oss << "      hint: \"DOF equation (index < nqsu=" << nqsu << "), not a constraint\"\n";
				} else {
					oss << "      hint: \"Unknown constraint type (not found in joints or parts)\"\n";
				}
			}
			oss << "---END:INCONSISTENT_CONSTRAINTS---\n";

			// Re-throw with detailed message
			throw InconsistentConstraintsError(oss.str());
		}
		catch (const SingularMatrixError& ex) {
			auto redundantEqnNos = ex.getRedundantEqnNos();
			auto rhsAtDetection = ex.getRhsValues();
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->removeRedundantConstraints(redundantEqnNos); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->constraintsReport(); });
			system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(qsuOld); });
			// Store equation numbers and RHS for post-convergence verification
			if (rhsAtDetection) {
				if (!removedEqnNos) {
					removedEqnNos = std::make_shared<std::vector<size_t>>();
					removedRhsAtDetection = std::make_shared<std::vector<double>>();
				}
				for (size_t i = 0; i < redundantEqnNos->size(); i++) {
					removedEqnNos->push_back(redundantEqnNos->at(i));
					removedRhsAtDetection->push_back(rhsAtDetection->at(i));
				}
			}
		}
		catch (SimulationStoppingError& ex) {
			// Check if this is a convergence failure
			std::string msg(ex.what());
			if (msg.find("iterNo > iterMax") != std::string::npos) {
				// Convergence failed - build diagnostic with joint information
				system->logString("---BEGIN:CONVERGENCE_FAILURE---");
				system->logString("Assembly solver failed to converge after maximum iterations");

				// Reactivate redundant constraints
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					item->reactivateRedundantConstraints();
				});

				// Collect all violated equation numbers and values
				auto violatedEqnNos = std::make_shared<FullColumn<size_t>>();
				auto violationValues = std::make_shared<std::vector<double>>();

				// Evaluate joint constraints
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					auto joint = std::dynamic_pointer_cast<Joint>(item);
					if (joint) {
						joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
							con->calcPostDynCorrectorIteration();
							double violation = std::abs(con->aG);
							if (violation > 1.0e-6) {
								violatedEqnNos->push_back(con->iG);
								violationValues->push_back(violation);
							}
						});
					}
				});

				// Evaluate part constraints
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					auto part = std::dynamic_pointer_cast<Part>(item);
					if (part && part->partFrame) {
						auto aGeu = part->partFrame->aGeu;
						if (aGeu) {
							aGeu->calcPostDynCorrectorIteration();
							double violation = std::abs(aGeu->aG);
							if (violation > 1.0e-6) {
								violatedEqnNos->push_back(aGeu->iG);
								violationValues->push_back(violation);
							}
						}
						if (part->partFrame->aGabs) {
							for (auto& aGab : *(part->partFrame->aGabs)) {
								aGab->calcPostDynCorrectorIteration();
								double violation = std::abs(aGab->aG);
								if (violation > 1.0e-6) {
									violatedEqnNos->push_back(aGab->iG);
									violationValues->push_back(violation);
								}
							}
						}
					}
				});

				// Build diagnostic using joint infrastructure
				std::ostringstream oss;
				oss << std::fixed << std::setprecision(6);
				oss << "  total_iterations: 101\n";
				oss << "  total_violations: " << violatedEqnNos->size() << "\n";

				double totalViolation = 0.0;
				for (auto v : *violationValues) {
					totalViolation += v;
				}
				oss << "  total_violation_magnitude: " << totalViolation << "\n";
				oss << "  joints:\n";

				// Get joint diagnostics (includes part names, LCS info, positions)
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					auto joint = std::dynamic_pointer_cast<Joint>(item);
					if (joint) {
						auto jointDiag = joint->getJointDiagnostic(violatedEqnNos, violationValues);
						if (!jointDiag.inconsistentConstraints.empty()) {
							oss << "    - name: \"" << jointDiag.name << "\"\n";
							oss << "      type: \"" << jointDiag.type << "\"\n";
							oss << "      part_i: \"" << jointDiag.partIName << "\"\n";
							oss << "      part_j: \"" << jointDiag.partJName << "\"\n";
							oss << "      lcs_i:\n";
							oss << "        name: \"" << jointDiag.lcsI.name << "\"\n";
							oss << "        world_position: [" << jointDiag.lcsI.worldPosition[0] << ", "
								<< jointDiag.lcsI.worldPosition[1] << ", " << jointDiag.lcsI.worldPosition[2] << "]\n";
							oss << "      lcs_j:\n";
							oss << "        name: \"" << jointDiag.lcsJ.name << "\"\n";
							oss << "        world_position: [" << jointDiag.lcsJ.worldPosition[0] << ", "
								<< jointDiag.lcsJ.worldPosition[1] << ", " << jointDiag.lcsJ.worldPosition[2] << "]\n";
							oss << "      violated_constraints:\n";
							for (const auto& conDiag : jointDiag.inconsistentConstraints) {
								oss << "        - type: \"" << conDiag.type << "\"\n";
								oss << "          violation: " << conDiag.violation << "\n";
								oss << "          equation: " << conDiag.equationNumber << "\n";
							}
						}
					}
				});

				// Add part constraints that don't belong to joints
				oss << "  part_constraints:\n";
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					auto part = std::dynamic_pointer_cast<Part>(item);
					if (part && part->partFrame) {
						auto aGeu = part->partFrame->aGeu;
						if (aGeu && violatedEqnNos) {
							auto it = std::find(violatedEqnNos->begin(), violatedEqnNos->end(), aGeu->iG);
							if (it != violatedEqnNos->end()) {
								size_t idx = std::distance(violatedEqnNos->begin(), it);
								double violation = (idx < violationValues->size()) ? violationValues->at(idx) : aGeu->aG;
								oss << "    - part: \"" << part->name << "\"\n";
								oss << "      constraint: \"" << aGeu->constraintSpec() << "\"\n";
								oss << "      violation: " << violation << "\n";
								oss << "      equation: " << aGeu->iG << "\n";
							}
						}
					}
				});

				system->logString(oss.str());
				system->logString("---END:CONVERGENCE_FAILURE---");

				// Re-throw to propagate error
				throw;
			} else {
				// Other SimulationStoppingError - re-throw
				throw;
			}
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

void PosICNewtonRaphson::verifyRemovedConstraintsAtConvergence()
{
	// If no constraints were removed, nothing to verify
	if (!removedEqnNos || removedEqnNos->empty()) {
		return;
	}

	// Get the current constraint residuals at the converged state
	// by asking each RedundantConstraint to compute its wrapped constraint's aG
	std::vector<std::pair<size_t, double>> inconsistentConstraints;
	double consistencyTolerance = 1.0e-6;

	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
		auto joint = std::dynamic_pointer_cast<Joint>(item);
		if (joint) {
			joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
				if (con->isRedundant()) {
					auto redunCon = std::static_pointer_cast<RedundantConstraint>(con);
					auto wrappedCon = redunCon->constraint;

					// Use postPosICIteration() which updates all intermediate helper objects
					// (like aAijIeJe in DirectionCosineConstraint) before computing aG.
					// Just calling calcPostDynCorrectorIteration() doesn't update the helpers.
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
			// Check Euler constraint
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
			// Check absolute constraints
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

	// If we found any constraints that are truly inconsistent at convergence, throw error
	if (!inconsistentConstraints.empty()) {
		// Reactivate constraints first so the diagnostic code can access them
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->reactivateRedundantConstraints(); });

		// Build the list of inconsistent equation numbers and RHS values
		auto inconsistentEqnNos = std::make_shared<FullColumn<size_t>>();
		auto rhsValues = std::make_shared<std::vector<double>>();
		for (const auto& pair : inconsistentConstraints) {
			inconsistentEqnNos->push_back(pair.first);
			rhsValues->push_back(pair.second);
		}

		throw InconsistentConstraintsError(
			"Constraints are geometrically inconsistent (no solution exists)", inconsistentEqnNos, rhsValues);
	}
}
