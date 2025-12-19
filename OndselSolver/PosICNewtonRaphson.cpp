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
#include "System.h"
#include "ExternalSystem.h"
#include "Part.h"
#include "PartFrame.h"
#include "Constraint.h"
#include "RedundantConstraint.h"
#include "CREATE.h"
#include "GESpMatParPvPrecise.h"
#include "GESpMatFullPvPosIC.h"
#include "Joint.h"
#include "EndFrameqc.h"
#include "MarkerFrame.h"

using namespace MbD;

// Helper: Hamilton product of two quaternions q1 × q2
// Quaternion format: (x, y, z, w) where w is scalar
static std::array<double, 4> quaternionMultiply(const std::array<double, 4>& q1, const std::array<double, 4>& q2) {
	double x1 = q1[0], y1 = q1[1], z1 = q1[2], w1 = q1[3];
	double x2 = q2[0], y2 = q2[1], z2 = q2[2], w2 = q2[3];
	return {
		w1*x2 + x1*w2 + y1*z2 - z1*y2,  // x
		w1*y2 - x1*z2 + y1*w2 + z1*x2,  // y
		w1*z2 + x1*y2 - y1*x2 + z1*w2,  // z
		w1*w2 - x1*x2 - y1*y2 - z1*z2   // w
	};
}

// Helper: Conjugate of quaternion (x, y, z, w) -> (-x, -y, -z, w)
static std::array<double, 4> quaternionConjugate(const std::array<double, 4>& q) {
	return {-q[0], -q[1], -q[2], q[3]};
}

void PosICNewtonRaphson::run()
{
	// Clear any previous tracking
	removedEqnNos = nullptr;
	removedRhsAtDetection = nullptr;
	protectedConstraints.clear();
	allRemovedConstraints.clear();
	correctedPartsFor180.clear();

	try {  // OUTER try - catches ALL InconsistentConstraintsError for YAML processing
		while (true) {
			try {
				//VectorNewtonRaphson::run();   //Inline to help debugging
				preRun();
				initializeLocally();
				initializeGlobally();
				iterate();
				postRun();
				// After successful convergence, verify removed constraints are satisfied
				// Returns true if essential constraints were found and protected for retry
				if (verifyRemovedConstraintsAtConvergence()) {
					continue;  // Retry the solve with protected constraints
				}
				break;
			}
			catch (const SingularMatrixError& ex) {
				auto redundantEqnNos = ex.getRedundantEqnNos();
				auto rhsAtDetection = ex.getRhsValues();

				// ============================================================
				// FIRST: Check for 180° configuration singularity
				// If detected, apply correction and retry instead of removing constraint
				// ============================================================
				bool applied180Correction = false;
				const double SINGULARITY_W_TOLERANCE = 0.1;  // |w| < 0.1 means angle > ~168°
				const double RESIDUAL_TOLERANCE = 0.1;       // Large residual = config singularity

				for (size_t eqnNo : *redundantEqnNos) {
					if (applied180Correction) break;

					// Find the joint containing this constraint and check for 180° singularity
					system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
						if (applied180Correction) return;

						auto joint = std::dynamic_pointer_cast<Joint>(item);
						if (!joint) return;

						// Check if this joint contains the singular constraint
						bool hasConstraint = false;
						double constraintResidual = 0.0;
						joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
							if (con->iG == eqnNo) {
								hasConstraint = true;
								// Get actual residual (aG is updated during iteration)
								constraintResidual = std::abs(con->aG);
							}
						});

						if (!hasConstraint) return;

						// Only proceed if residual is large (configuration singularity, not redundancy)
						if (constraintResidual < RESIDUAL_TOLERANCE) return;

						// Get end frames as EndFrameqc to access world quaternions
						auto frmIqc = std::dynamic_pointer_cast<EndFrameqc>(joint->frmI);
						auto frmJqc = std::dynamic_pointer_cast<EndFrameqc>(joint->frmJ);
						if (!frmIqc || !frmJqc) return;

						// Get world quaternions (accounts for marker frames)
						auto qI_ptr = frmIqc->qEO();
						auto qJ_ptr = frmJqc->qEO();

						std::array<double, 4> qI = {qI_ptr->at(0), qI_ptr->at(1), qI_ptr->at(2), qI_ptr->at(3)};
						std::array<double, 4> qJ = {qJ_ptr->at(0), qJ_ptr->at(1), qJ_ptr->at(2), qJ_ptr->at(3)};

						// Compute relative quaternion: q_rel = conj(qI) × qJ
						auto qI_conj = quaternionConjugate(qI);
						auto qRel = quaternionMultiply(qI_conj, qJ);

						// Check if near 180° rotation (pure imaginary quaternion: |w| ≈ 0)
						if (std::abs(qRel[3]) >= SINGULARITY_W_TOLERANCE) return;

						// This is a 180° configuration singularity!
						// Find dominant axis and construct correction quaternion
						double absX = std::abs(qRel[0]);
						double absY = std::abs(qRel[1]);
						double absZ = std::abs(qRel[2]);

						std::array<double, 4> qCorr = {0, 0, 0, 0};
						std::string axisName;
						if (absX >= absY && absX >= absZ) {
							qCorr[0] = (qRel[0] >= 0) ? 1.0 : -1.0;
							axisName = "x";
						} else if (absY >= absZ) {
							qCorr[1] = (qRel[1] >= 0) ? 1.0 : -1.0;
							axisName = "y";
						} else {
							qCorr[2] = (qRel[2] >= 0) ? 1.0 : -1.0;
							axisName = "z";
						}

						// Get part J's quaternion and apply correction
						auto markerJ = frmJqc->getMarkerFrame();
						if (!markerJ || !markerJ->partFrame) return;

						auto partFrame = markerJ->partFrame;

						// Check if this part was already corrected (prevents infinite loop with shared parts)
						if (correctedPartsFor180.count(partFrame) > 0) {
							// This part was already corrected by another joint - skip
							return;
						}

						auto qE_old = partFrame->qE;

						// Get marker J's local quaternion (rotation from part frame to marker frame)
						// qEpm is the quaternion form of aApm (constant, computed from marker placement)
						auto qMarkerJ = markerJ->qEpm;
						std::array<double, 4> qM = {qMarkerJ->at(0), qMarkerJ->at(1), qMarkerJ->at(2), qMarkerJ->at(3)};
						auto qM_conj = quaternionConjugate(qM);

						// Conjugate correction by marker: qCorr_adj = qMarker × qCorr × conj(qMarker)
						// This transforms the end-frame correction to a part-frame correction
						// Math: qEndFrame' = qPart' × qMarker = qPart × qMarker × qCorr
						//       So: qPart' = qPart × qMarker × qCorr × conj(qMarker)
						auto qTemp = quaternionMultiply(qM, qCorr);
						auto qCorr_adj = quaternionMultiply(qTemp, qM_conj);

						// Compute new quaternion: qE_new = qE_old × qCorr_adj
						std::array<double, 4> qOld = {qE_old->at(0), qE_old->at(1), qE_old->at(2), qE_old->at(3)};
						auto qNew = quaternionMultiply(qOld, qCorr_adj);

						// Normalize for numerical stability
						double norm = std::sqrt(qNew[0]*qNew[0] + qNew[1]*qNew[1] + qNew[2]*qNew[2] + qNew[3]*qNew[3]);
						if (norm > 1e-10) {
							qNew[0] /= norm;
							qNew[1] /= norm;
							qNew[2] /= norm;
							qNew[3] /= norm;
						}

						// Apply correction to part's quaternion
						qE_old->at(0) = qNew[0];
						qE_old->at(1) = qNew[1];
						qE_old->at(2) = qNew[2];
						qE_old->at(3) = qNew[3];

						// Log the correction
						std::ostringstream msg;
						msg << "MbD: Detected 180° configuration singularity in joint '" << joint->name
						    << "' (relative_w=" << std::abs(qRel[3]) << ", residual=" << constraintResidual
						    << "). Applied 180° correction about " << axisName << "-axis to part '"
						    << partFrame->part->name << "'.";
						system->logString(msg.str());

						// Mark part as corrected to prevent re-correction by other joints
						correctedPartsFor180.insert(partFrame);

						applied180Correction = true;
					});
				}

				if (applied180Correction) {
					// Retry iteration from corrected state (don't reset to qsuOld)
					continue;
				}

				// ============================================================
				// No 180° singularity detected - proceed with normal redundancy handling
				// ============================================================

				// Helper to check if a constraint (by its current iG) is protected
				// We need to find the constraint object for a given equation number
				auto isEqnProtected = [&](size_t eqnNo) -> bool {
					bool found = false;
					system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
						if (found) return;
						auto joint = std::dynamic_pointer_cast<Joint>(item);
						if (joint) {
							joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
								if (con->iG == eqnNo && protectedConstraints.count(con.get()) > 0) {
									found = true;
								}
							});
						}
						auto part = std::dynamic_pointer_cast<Part>(item);
						if (part && part->partFrame) {
							if (part->partFrame->aGeu && part->partFrame->aGeu->iG == eqnNo &&
								protectedConstraints.count(part->partFrame->aGeu.get()) > 0) {
								found = true;
							}
							if (part->partFrame->aGabs) {
								for (auto& aGab : *(part->partFrame->aGabs)) {
									if (aGab->iG == eqnNo && protectedConstraints.count(aGab.get()) > 0) {
										found = true;
									}
								}
							}
						}
					});
					return found;
				};

				// Filter out protected constraints - these were learned to be essential in previous iterations
				auto toRemove = std::make_shared<FullColumn<size_t>>();
				auto toRemoveRhs = std::make_shared<std::vector<double>>();

				for (size_t i = 0; i < redundantEqnNos->size(); i++) {
					size_t eqnNo = redundantEqnNos->at(i);
					if (!isEqnProtected(eqnNo)) {
						toRemove->push_back(eqnNo);
						if (rhsAtDetection && i < rhsAtDetection->size()) {
							toRemoveRhs->push_back(rhsAtDetection->at(i));
						}
					}
				}

				// If all singular constraints are protected, we're stuck - can't converge, can't remove
				if (toRemove->empty()) {
					// Build list of protected constraint equation numbers
					auto eqnNosFullCol = std::make_shared<FullColumn<size_t>>();
					for (Constraint* con : protectedConstraints) {
						eqnNosFullCol->push_back(con->iG);
					}

					// Outer catch will: reset to qsuOld, compute residuals, build YAML,
					// restore bestEffortState, then rethrow
					throw InconsistentConstraintsError(
						"Cannot solve: all singular constraints are protected (geometrically inconsistent)",
						eqnNosFullCol, rhsAtDetection);
				}

				// Track constraint pointers before they're wrapped (for final violation reporting)
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
					auto joint = std::dynamic_pointer_cast<Joint>(item);
					if (joint) {
						joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
							if (std::find(toRemove->begin(), toRemove->end(), con->iG) != toRemove->end()) {
								allRemovedConstraints.insert(con.get());
							}
						});
					}
					auto part = std::dynamic_pointer_cast<Part>(item);
					if (part && part->partFrame) {
						if (part->partFrame->aGeu) {
							if (std::find(toRemove->begin(), toRemove->end(), part->partFrame->aGeu->iG) != toRemove->end()) {
								allRemovedConstraints.insert(part->partFrame->aGeu.get());
							}
						}
						if (part->partFrame->aGabs) {
							for (auto& aGab : *(part->partFrame->aGabs)) {
								if (std::find(toRemove->begin(), toRemove->end(), aGab->iG) != toRemove->end()) {
									allRemovedConstraints.insert(aGab.get());
								}
							}
						}
					}
				});

				// Remove only non-protected constraints
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->removeRedundantConstraints(toRemove); });
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->constraintsReport(); });
				system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(qsuOld); });

				// Store equation numbers and RHS for post-convergence verification
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
							auto jointDiag = joint->getJointDiagnostic(violatedEqnNos);
							if (!jointDiag.inconsistentConstraints.empty()) {
								// Helper to convert quaternion [x,y,z,w] to Euler angles [roll,pitch,yaw] in degrees
								auto quatToEuler = [](const std::array<double, 4>& q) -> std::array<double, 3> {
									double x = q[0], y = q[1], z = q[2], w = q[3];
									double roll = std::atan2(2.0*(w*x + y*z), 1.0 - 2.0*(x*x + y*y)) * 180.0 / M_PI;
									double sinp = 2.0*(w*y - z*x);
									double pitch = (std::abs(sinp) >= 1.0) ? std::copysign(90.0, sinp) : std::asin(sinp) * 180.0 / M_PI;
									double yaw = std::atan2(2.0*(w*z + x*y), 1.0 - 2.0*(y*y + z*z)) * 180.0 / M_PI;
									return {roll, pitch, yaw};
								};
								auto eulerI = quatToEuler(jointDiag.lcsI.worldQuaternion);
								auto eulerJ = quatToEuler(jointDiag.lcsJ.worldQuaternion);

								oss << "    - type: \"" << jointDiag.type << "\"\n";
								oss << "      part_i:\n";
								oss << "        name: \"" << jointDiag.partIName << "\"\n";
								oss << "        lcs_position_on_part: [" << jointDiag.lcsI.positionOnPart[0] << ", "
									<< jointDiag.lcsI.positionOnPart[1] << ", " << jointDiag.lcsI.positionOnPart[2] << "]\n";
								oss << "        lcs_world_position: [" << jointDiag.lcsI.worldPosition[0] << ", "
									<< jointDiag.lcsI.worldPosition[1] << ", " << jointDiag.lcsI.worldPosition[2] << "]\n";
								oss << "        lcs_world_quaternion: [" << jointDiag.lcsI.worldQuaternion[0] << ", "
									<< jointDiag.lcsI.worldQuaternion[1] << ", " << jointDiag.lcsI.worldQuaternion[2] << ", "
									<< jointDiag.lcsI.worldQuaternion[3] << "]\n";
								oss << "        lcs_world_orientation: [" << eulerI[0] << ", " << eulerI[1] << ", " << eulerI[2] << "]\n";
								oss << "      part_j:\n";
								oss << "        name: \"" << jointDiag.partJName << "\"\n";
								oss << "        lcs_position_on_part: [" << jointDiag.lcsJ.positionOnPart[0] << ", "
									<< jointDiag.lcsJ.positionOnPart[1] << ", " << jointDiag.lcsJ.positionOnPart[2] << "]\n";
								oss << "        lcs_world_position: [" << jointDiag.lcsJ.worldPosition[0] << ", "
									<< jointDiag.lcsJ.worldPosition[1] << ", " << jointDiag.lcsJ.worldPosition[2] << "]\n";
								oss << "        lcs_world_quaternion: [" << jointDiag.lcsJ.worldQuaternion[0] << ", "
									<< jointDiag.lcsJ.worldQuaternion[1] << ", " << jointDiag.lcsJ.worldQuaternion[2] << ", "
									<< jointDiag.lcsJ.worldQuaternion[3] << "]\n";
								oss << "        lcs_world_orientation: [" << eulerJ[0] << ", " << eulerJ[1] << ", " << eulerJ[2] << "]\n";
								oss << "      relative_angle_degrees: " << jointDiag.relativeAngleDegrees << "\n";
								oss << "      relative_quaternion: [" << jointDiag.relativeQuaternion[0] << ", "
									<< jointDiag.relativeQuaternion[1] << ", " << jointDiag.relativeQuaternion[2] << ", "
									<< jointDiag.relativeQuaternion[3] << "]\n";
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
									// Use aGeu->aG directly - calcPostDynCorrectorIteration() was called above
									oss << "    - part: \"" << part->name << "\"\n";
									oss << "      constraint: \"" << aGeu->constraintSpec() << "\"\n";
									oss << "      violation: " << aGeu->aG << "\n";
									oss << "      equation: " << aGeu->iG << "\n";
								}
							}
						}
					});

					system->logString(oss.str());
					system->logString("---END:CONVERGENCE_FAILURE---");

					// Restore best state seen during iteration (lowest yNorm)
					// This gives user the closest-to-solution configuration, not the diverged final state
					if (bestEffortState) {
						system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
							item->setqsu(bestEffortState);
						});
					}

					// Update parts with best-effort positions for visualization
					system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
						auto part = std::dynamic_pointer_cast<Part>(item);
						if (part) {
							part->postPosICIteration();
						}
					});

					// Push best-effort state to FreeCAD
					system->system->externalSystem->updateFromMbD();
					return;  // Return normally - don't re-throw
				} else {
					// Other SimulationStoppingError - re-throw
					throw;
				}
			}
		}
	}
	catch (InconsistentConstraintsError& ex) {
		// ALL InconsistentConstraintsError flow here - build YAML diagnostic
		auto inconsistentEqnNos = ex.getInconsistentEqnNos();
		// NOTE: We ignore ex.getRhsValues() - those are from Gaussian elimination, not actual residuals

		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->reactivateRedundantConstraints(); });

		// Use best-effort state (saved in verifyRemovedConstraintsAtConvergence before retry)
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsu(bestEffortState); });

		// Update derived quantities (rOeO, etc.) on Parts and their marker frames
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			auto part = std::dynamic_pointer_cast<Part>(item);
			if (part) {
				part->postPosICIteration();
			}
		});

		// Compute ACTUAL constraint residuals at current state by calling postPosICIteration()
		// This updates each constraint's aG value to the true residual
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			auto joint = std::dynamic_pointer_cast<Joint>(item);
			if (joint) {
				joint->constraintsDo([&](std::shared_ptr<Constraint> con) {
					con->postPosICIteration();  // Updates con->aG to actual residual
				});
			}
		});
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			auto part = std::dynamic_pointer_cast<Part>(item);
			if (part && part->partFrame) {
				if (part->partFrame->aGeu) {
					part->partFrame->aGeu->postPosICIteration();
				}
				if (part->partFrame->aGabs) {
					for (auto& aGab : *(part->partFrame->aGabs)) {
						aGab->postPosICIteration();
					}
				}
			}
		});

			// Build diagnostic information
			auto diagnostic = std::make_shared<InconsistencyDiagnostic>();
			std::map<std::string, int> partOccurrences;
			std::set<size_t> matchedEqnNos;  // Track which equation numbers were matched

			// Collect joint diagnostics
		// getJointDiagnostic uses con->aG (actual residual computed above by postPosICIteration)
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			auto joint = std::dynamic_pointer_cast<Joint>(item);
			if (joint) {
				auto jointDiag = joint->getJointDiagnostic(inconsistentEqnNos);
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
		// Use aG directly since we called postPosICIteration() above
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
						pcDiag.violation = aGeu->aG;  // Use actual residual
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
							pcDiag.violation = aGab->aG;  // Use actual residual
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
					// We don't have residual for unmatched equations, use 0.0
					unmatchedEqns.push_back({eqnNo, 0.0});
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
			// Helper to convert quaternion [x,y,z,w] to Euler angles [roll,pitch,yaw] in degrees
			auto quatToEuler = [](const std::array<double, 4>& q) -> std::array<double, 3> {
				double x = q[0], y = q[1], z = q[2], w = q[3];
				double roll = std::atan2(2.0*(w*x + y*z), 1.0 - 2.0*(x*x + y*y)) * 180.0 / M_PI;
				double sinp = 2.0*(w*y - z*x);
				double pitch = (std::abs(sinp) >= 1.0) ? std::copysign(90.0, sinp) : std::asin(sinp) * 180.0 / M_PI;
				double yaw = std::atan2(2.0*(w*z + x*y), 1.0 - 2.0*(y*y + z*z)) * 180.0 / M_PI;
				return {roll, pitch, yaw};
			};
			for (const auto& jointDiag : diagnostic->joints) {
				auto eulerI = quatToEuler(jointDiag.lcsI.worldQuaternion);
				auto eulerJ = quatToEuler(jointDiag.lcsJ.worldQuaternion);

				oss << "    - type: \"" << jointDiag.type << "\"\n";
				oss << "      part_i:\n";
				oss << "        name: \"" << jointDiag.partIName << "\"\n";
				oss << "        lcs_position_on_part: [" << jointDiag.lcsI.positionOnPart[0] << ", "
					<< jointDiag.lcsI.positionOnPart[1] << ", " << jointDiag.lcsI.positionOnPart[2] << "]\n";
				oss << "        lcs_world_position: [" << jointDiag.lcsI.worldPosition[0] << ", "
					<< jointDiag.lcsI.worldPosition[1] << ", " << jointDiag.lcsI.worldPosition[2] << "]\n";
				oss << "        lcs_world_quaternion: [" << jointDiag.lcsI.worldQuaternion[0] << ", "
					<< jointDiag.lcsI.worldQuaternion[1] << ", " << jointDiag.lcsI.worldQuaternion[2] << ", "
					<< jointDiag.lcsI.worldQuaternion[3] << "]\n";
				oss << "        lcs_world_orientation: [" << eulerI[0] << ", " << eulerI[1] << ", " << eulerI[2] << "]\n";
				oss << "      part_j:\n";
				oss << "        name: \"" << jointDiag.partJName << "\"\n";
				oss << "        lcs_position_on_part: [" << jointDiag.lcsJ.positionOnPart[0] << ", "
					<< jointDiag.lcsJ.positionOnPart[1] << ", " << jointDiag.lcsJ.positionOnPart[2] << "]\n";
				oss << "        lcs_world_position: [" << jointDiag.lcsJ.worldPosition[0] << ", "
					<< jointDiag.lcsJ.worldPosition[1] << ", " << jointDiag.lcsJ.worldPosition[2] << "]\n";
				oss << "        lcs_world_quaternion: [" << jointDiag.lcsJ.worldQuaternion[0] << ", "
					<< jointDiag.lcsJ.worldQuaternion[1] << ", " << jointDiag.lcsJ.worldQuaternion[2] << ", "
					<< jointDiag.lcsJ.worldQuaternion[3] << "]\n";
				oss << "        lcs_world_orientation: [" << eulerJ[0] << ", " << eulerJ[1] << ", " << eulerJ[2] << "]\n";
				oss << "      relative_angle_degrees: " << jointDiag.relativeAngleDegrees << "\n";
				oss << "      relative_quaternion: [" << jointDiag.relativeQuaternion[0] << ", "
					<< jointDiag.relativeQuaternion[1] << ", " << jointDiag.relativeQuaternion[2] << ", "
					<< jointDiag.relativeQuaternion[3] << "]\n";
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

		// Log the diagnostic YAML
		system->logString(oss.str());

		// Push best-effort state to FreeCAD (already set at start of handler)
		system->system->externalSystem->updateFromMbD();
		return;
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
	// Track best state seen during iteration (lowest yNorm = best solution)
	// This is called after every iteration, right after yNorm is updated
	if (iterNo == 0) {
		bestYNorm = std::numeric_limits<double>::max();
		if (nqsu > 0) {
			bestEffortState = std::make_shared<FullColumn<double>>(nqsu);
		}
	}

	// Save state if this is the best we've seen (lower yNorm = better)
	if (bestEffortState && yNorm < bestYNorm) {
		bestYNorm = yNorm;
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			item->fillqsu(bestEffortState);
		});
	}

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

bool PosICNewtonRaphson::verifyRemovedConstraintsAtConvergence()
{
	// If no constraints were removed, nothing to verify
	if (!removedEqnNos || removedEqnNos->empty()) {
		return false;  // No retry needed
	}

	// Get the current constraint residuals at the converged state
	// by asking each RedundantConstraint to compute its wrapped constraint's aG
	// Store constraint pointers (not equation numbers) because iG changes on each retry
	std::vector<std::pair<Constraint*, double>> essentialConstraints;
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
					wrappedCon->postPosICIteration();
					double residual = wrappedCon->aG;
					if (std::abs(residual) > consistencyTolerance) {
						// This constraint has a non-zero residual at convergence
						// meaning it was essential and should NOT have been removed
						// Check for duplicates before adding
						auto it = std::find_if(essentialConstraints.begin(), essentialConstraints.end(),
							[&](const auto& p) { return p.first == wrappedCon.get(); });
						if (it == essentialConstraints.end()) {
							essentialConstraints.push_back({wrappedCon.get(), residual});
						}
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
					// Check for duplicates before adding
					auto it = std::find_if(essentialConstraints.begin(), essentialConstraints.end(),
						[&](const auto& p) { return p.first == wrappedCon.get(); });
					if (it == essentialConstraints.end()) {
						essentialConstraints.push_back({wrappedCon.get(), residual});
					}
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
							// Check for duplicates before adding
							auto it = std::find_if(essentialConstraints.begin(), essentialConstraints.end(),
								[&](const auto& p) { return p.first == wrappedCon.get(); });
							if (it == essentialConstraints.end()) {
								essentialConstraints.push_back({wrappedCon.get(), residual});
							}
						}
					}
				}
			}
		}
	});

	if (!essentialConstraints.empty()) {
		// These constraints were incorrectly removed - they are essential, not redundant
		// Protect them from future removal and retry the solve
		std::ostringstream protectedMsg;
		protectedMsg << "MbD: Protecting " << essentialConstraints.size()
		             << " essential constraint(s) from removal:";
		for (const auto& pair : essentialConstraints) {
			protectedConstraints.insert(pair.first);
			protectedMsg << " " << pair.first->iG;
		}
		system->logString(protectedMsg.str());

		// Save current converged state as "best effort" before resetting
		// This will be restored if we ultimately detect inconsistent constraints
		bestEffortState = std::make_shared<FullColumn<double>>(nqsu);
		system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
			item->fillqsu(bestEffortState);
		});

		// Reactivate ALL constraints and retry with protections in place
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

		return true;  // Signal to retry with protected constraints
	}

	return false;  // No retry needed - all removed constraints are truly redundant
}
