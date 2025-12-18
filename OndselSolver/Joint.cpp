/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <algorithm>
#include <memory>
#include <typeinfo>
#include <cmath>

#include "Joint.h"
#include "Constraint.h"
#include "EndFrameqc.h"
#include "EndFrameqct.h"
#include "CREATE.h"
#include "RedundantConstraint.h"
#include "MarkerFrame.h"
#include "ForceTorqueData.h"
#include "System.h"
#include "PartFrame.h"
#include "Part.h"

using namespace MbD;

Joint::Joint() {

}

Joint::Joint(const std::string& str) : ConstraintSet(str) {

}

void Joint::initializeLocally()
{
	auto frmIqc = std::dynamic_pointer_cast<EndFrameqc>(frmI);
	if (frmIqc) {
		if (frmIqc->endFrameqct) {
			frmI = frmIqc->endFrameqct;
		}
	}
	constraintsDo([](std::shared_ptr<Constraint> constraint) { constraint->initializeLocally(); });
}

FColDsptr MbD::Joint::aFIeJtIe()
{
	//"aFIeJtIe is joint force on end frame Ie expresses in Ie components."
	auto frmIqc = std::dynamic_pointer_cast<EndFrameqc>(frmI);
	return frmIqc->aAeO()->timesFullColumn(this->aFIeJtO());
}

FColDsptr MbD::Joint::aFIeJtO()
{
	//"aFIeJtO is joint force on end frame Ie expresses in O components."
	auto aFIeJtO = std::make_shared <FullColumn<double>>(3);
	constraintsDo([&](std::shared_ptr<Constraint> con) { con->addToJointForceI(aFIeJtO); });
	return aFIeJtO;
}

void Joint::fillRedundantConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> redunConstraints)
{
	constraintsDo([&](std::shared_ptr<Constraint> con) { con->fillRedundantConstraints(con, redunConstraints); });
}

void Joint::removeRedundantConstraints(std::shared_ptr<std::vector<size_t>> redundantEqnNos)
{
	for (size_t i = 0; i < constraints->size(); i++)
	{
		auto& constraint = constraints->at(i);
		if (std::find(redundantEqnNos->begin(), redundantEqnNos->end(), constraint->iG) != redundantEqnNos->end()) {
			auto redunCon = CREATE<RedundantConstraint>::With();
			redunCon->constraint = constraint;
			constraints->at(i) = redunCon;
		}
	}
}

void Joint::reactivateRedundantConstraints()
{
	for (size_t i = 0; i < constraints->size(); i++)
	{
		auto& con = constraints->at(i);
		if (con->isRedundant()) {
			constraints->at(i) = std::static_pointer_cast<RedundantConstraint>(con)->constraint;
		}
	}
}

void Joint::constraintsReport()
{
	auto redunCons = std::make_shared<std::vector<std::shared_ptr<Constraint>>>();
	constraintsDo([&](std::shared_ptr<Constraint> con) {
		if (con->isRedundant()) {
			redunCons->push_back(con);
		}
		});
	if (redunCons->size() > 0) {
		std::string str = "MbD: " + this->classname() + std::string(" ") + this->name + " has the following constraint(s) removed: ";
		this->logString(str);
		std::for_each(redunCons->begin(), redunCons->end(), [&](auto& con) {
			str = "MbD: " + std::string("    ") + con->constraintSpec();
			this->logString(str);
			});
	}
}

void Joint::inconsistentConstraintsReport(std::shared_ptr<std::vector<size_t>> inconsistentEqnNos)
{
	auto inconsistentCons = std::make_shared<std::vector<std::shared_ptr<Constraint>>>();
	constraintsDo([&](std::shared_ptr<Constraint> con) {
		if (std::find(inconsistentEqnNos->begin(), inconsistentEqnNos->end(), con->iG) != inconsistentEqnNos->end()) {
			inconsistentCons->push_back(con);
		}
		});
	if (inconsistentCons->size() > 0) {
		std::string str = "MbD: " + this->classname() + std::string(" ") + this->name + " has the following inconsistent constraint(s): ";
		this->logString(str);
		std::for_each(inconsistentCons->begin(), inconsistentCons->end(), [&](auto& con) {
			str = "MbD: " + std::string("    ") + con->constraintSpec();
			this->logString(str);
			});
	}
}

JointDiagnostic Joint::getJointDiagnostic(std::shared_ptr<std::vector<size_t>> inconsistentEqnNos)
{
	JointDiagnostic diag;
	diag.name = this->name;
	diag.type = this->classname();

	// Get part and LCS info from frmI
	auto markerI = frmI->getMarkerFrame();
	auto partFrameI = markerI->getPartFrame();
	auto partI = partFrameI->getPart();
	diag.partIName = partI ? partI->name : "";
	diag.lcsI.name = markerI->name;
	diag.lcsI.positionOnPart = { markerI->rpmp->at(0), markerI->rpmp->at(1), markerI->rpmp->at(2) };
	diag.lcsI.worldPosition = { frmI->rOeO->at(0), frmI->rOeO->at(1), frmI->rOeO->at(2) };

	// Get part and LCS info from frmJ
	auto markerJ = frmJ->getMarkerFrame();
	auto partFrameJ = markerJ->getPartFrame();
	auto partJ = partFrameJ->getPart();
	diag.partJName = partJ ? partJ->name : "";
	diag.lcsJ.name = markerJ->name;
	diag.lcsJ.positionOnPart = { markerJ->rpmp->at(0), markerJ->rpmp->at(1), markerJ->rpmp->at(2) };
	diag.lcsJ.worldPosition = { frmJ->rOeO->at(0), frmJ->rOeO->at(1), frmJ->rOeO->at(2) };

	// Get world orientations
	auto efrmI = std::dynamic_pointer_cast<EndFrameqc>(frmI);
	auto efrmJ = std::dynamic_pointer_cast<EndFrameqc>(frmJ);
	if (efrmI && efrmJ) {
		auto qI = efrmI->qEO();
		auto qJ = efrmJ->qEO();

		diag.lcsI.worldQuaternion = {qI->at(0), qI->at(1), qI->at(2), qI->at(3)};
		diag.lcsJ.worldQuaternion = {qJ->at(0), qJ->at(1), qJ->at(2), qJ->at(3)};

		// Compute relative quaternion: conj(qI) * qJ
		// conj(qI) = (-qI0, -qI1, -qI2, qI3)
		double relW = qI->at(3)*qJ->at(3) + qI->at(0)*qJ->at(0)
		            + qI->at(1)*qJ->at(1) + qI->at(2)*qJ->at(2);
		double relX = qI->at(3)*qJ->at(0) - qI->at(0)*qJ->at(3)
		            - qI->at(1)*qJ->at(2) + qI->at(2)*qJ->at(1);
		double relY = qI->at(3)*qJ->at(1) + qI->at(0)*qJ->at(2)
		            - qI->at(1)*qJ->at(3) - qI->at(2)*qJ->at(0);
		double relZ = qI->at(3)*qJ->at(2) - qI->at(0)*qJ->at(1)
		            + qI->at(1)*qJ->at(0) - qI->at(2)*qJ->at(3);

		diag.relativeQuaternion = {relX, relY, relZ, relW};

		// Angle = 2 * acos(|w|) in degrees
		double clampedW = std::fabs(relW);
		if (clampedW > 1.0) clampedW = 1.0;
		diag.relativeAngleDegrees = 2.0 * std::acos(clampedW) * 180.0 / M_PI;
	}

	// Get inconsistent constraints for this joint
	// Always use con->aG - caller must call postPosICIteration() first to compute actual residuals
	constraintsDo([&](std::shared_ptr<Constraint> con) {
		auto it = std::find(inconsistentEqnNos->begin(), inconsistentEqnNos->end(), con->iG);
		if (it != inconsistentEqnNos->end()) {
			ConstraintDiagnostic conDiag;
			conDiag.equationNumber = con->iG;
			conDiag.type = con->constraintSpec();
			conDiag.violation = con->aG;
			diag.inconsistentConstraints.push_back(conDiag);
		}
		});

	return diag;
}

std::shared_ptr<StateData> Joint::stateData()
{
	//"
	//MbD returns aFIeO and aTIeO.
	//GEO needs aFImO and aTImO.
	//For Motion rImIeO is not zero and is changing.
	//aFImO : = aFIeO.
	//aTImO : = aTIeO + (rImIeO cross : aFIeO).
	//"

	auto answer = std::make_shared<ForceTorqueData>();
	auto aFIeO = this->aFX();
	auto aTIeO = this->aTX();
	auto rImIeO = this->frmI->rmeO();
	answer->aFIO = aFIeO;
	answer->aTIO = aTIeO->plusFullColumn(rImIeO->cross(aFIeO));
	return answer;
}

FColDsptr Joint::aFX()
{
	return this->jointForceI();
}

FColDsptr MbD::Joint::aTIeJtIe()
{
	//"aTIeJtIe is torque on part containing end frame Ie expressed in Ie components."
	return frmI->aAeO()->timesFullColumn(this->aTIeJtO());
}

FColDsptr MbD::Joint::aTIeJtO()
{
	//"aTIeJtO is torque on part containing end frame Ie expressed in O components."
	auto aTIeJtO = std::make_shared <FullColumn<double>>(3);
	constraintsDo([&](std::shared_ptr<Constraint> con) { con->addToJointTorqueI(aTIeJtO); });
	return aTIeJtO;
}

FColDsptr Joint::jointForceI()
{
	//"jointForceI is force on MbD marker I."
	auto jointForce = std::make_shared <FullColumn<double>>(3);
	constraintsDo([&](std::shared_ptr<Constraint> con) { con->addToJointForceI(jointForce); });
	return jointForce;
}

FColDsptr Joint::aTX()
{
	return this->jointTorqueI();
}

FColDsptr Joint::jointTorqueI()
{
	//"jointTorqueI is torque on MbD marker I."
	auto jointTorque = std::make_shared <FullColumn<double>>(3);
	constraintsDo([&](std::shared_ptr<Constraint> con) { con->addToJointTorqueI(jointTorque); });
	return jointTorque;
}
